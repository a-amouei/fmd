/*
  integrators.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2022 Arham Amouye Foumani

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <tgmath.h>
#include "fmd-private.h"
#include "matter.h"
#include "misc.h"
#include "turi.h"
#include "forces.h"
#include "timer.h"
#include "general.h"
#include "ttm.h"

#define WHAT_IF_THE_PARTICLE_HAS_LEFT_THE_BOX(md, c, i, d, x)                  \
    if ( ((md)->ns[(d)] == 1) && (((x) < 0.0) || ((x) >= (md)->l[(d)])) )      \
    {                                                                          \
        if (!(md)->PBC[(d)])                                                   \
        {                                                                      \
            _fmd_cell_remove_atom((md), (c), (i));                             \
            (md)->Subdomain.NumberOfParticles--;                               \
            (i)--;                                                             \
            break;                                                             \
        }                                                                      \
        else                                                                   \
            (x) += ((x) < 0.0 ? (md)->l[(d)] : -(md)->l[(d)]);                 \
    }                                                                          \
    do {} while (0)

fmd_real_t fmd_dync_getTimeStep(fmd_t *md)
{
    return md->timestep;
}

fmd_real_t fmd_dync_getTime(fmd_t *md)
{
    return md->time;
}

static void VelocityVerlet_startStep(fmd_t *md, fmd_bool_t UseThermostat)
{
    fmd_ituple_t ic;
    fmd_real_t VelocityScale, mass;
    cell_t *c;
    unsigned i;

    if (UseThermostat)
        VelocityScale = sqrt(1 + md->timestep / md->BerendsenThermostatParam *
                        (md->DesiredTemperature / md->GroupTemperature - 1));

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            for (int d=0; d<DIM; d++)
            {
                if (UseThermostat) VEL(c, i, d) *= VelocityScale;

                VEL(c, i, d) += md->timestep * 0.5 / mass * FRC(c, i, d);
                POS(c, i, d) += md->timestep * VEL(c, i, d);

                WHAT_IF_THE_PARTICLE_HAS_LEFT_THE_BOX(md, c, i, d, POS(c, i, d));
            }
        }

    /* check if the particles must be placed in other cells or subdomains */
    _fmd_refreshGrid(md);
}

/* updates GroupMomentum & GroupParticlesNum & GroupKineticEnergy & GroupTemperature */
void updateGroupProperties(fmd_t *md, fmd_rtuple_t MomentumSum, int ParticlesNum, fmd_real_t m_vSqd_Sum)
{
    fmd_real_t m_vSqd_SumSum;

    MPI_Reduce(MomentumSum, md->GroupMomentum, DIM, FMD_MPI_REAL, MPI_SUM, RANK0, md->MD_comm);

    MPI_Allreduce(&ParticlesNum, &md->GroupParticlesNum, 1, MPI_INT, MPI_SUM, md->MD_comm);

    MPI_Allreduce(&m_vSqd_Sum, &m_vSqd_SumSum, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
    md->GroupKineticEnergy = 0.5 * m_vSqd_SumSum;

    md->GroupTemperature = m_vSqd_SumSum / (3.0 * md->GroupParticlesNum * K_BOLTZMANN);
}

static void VelocityVerlet_finishStep(fmd_t *md)
{
    fmd_ituple_t ic;
    fmd_real_t m_vSqd_Sum = 0;
    fmd_real_t mass;
    int ParticlesNum = 0;
    fmd_rtuple_t MomentumSum = {0., 0., 0.};
    cell_t *c;
    unsigned i;

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            ParticlesNum++;

            mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            for (int d=0; d<DIM; d++)
            {
                VEL(c, i, d) += md->timestep * 0.5 / mass * FRC(c, i, d);
                MomentumSum[d] += mass * VEL(c, i, d);
            }

            m_vSqd_Sum += mass * ( sqrr(VEL(c, i, 0)) +
                                   sqrr(VEL(c, i, 1)) +
                                   sqrr(VEL(c, i, 2)) );
        }

    /* update GroupMomentum & GroupParticlesNum & GroupKineticEnergy & GroupTemperature */
    updateGroupProperties(md, MomentumSum, ParticlesNum, m_vSqd_Sum);
}

/* implementation of semi-implicit Euler method */
static void SymplecticEuler_takeOneStep(fmd_t *md)
{
    fmd_ituple_t ic;
    int d;
    fmd_real_t mass;
    cell_t *c;
    unsigned i;
    int ParticlesNum = 0;
    fmd_rtuple_t MomentumSum = {0., 0., 0.};
    fmd_real_t m_vSqd_Sum = 0;

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            for (d=0; d<DIM; d++)
            {
                VEL(c, i, d) += md->timestep * FRC(c, i, d) / mass;
                POS(c, i, d) += md->timestep * VEL(c, i, d);

                WHAT_IF_THE_PARTICLE_HAS_LEFT_THE_BOX(md, c, i, d, POS(c, i, d));
            }

            if (d == DIM)  /* if particle is not removed */
            {
                ParticlesNum++;

                fmd_real_t temp = 0.0;

                for (d=0; d<DIM; d++)
                {
                    temp += sqrr(VEL(c, i, d));
                    MomentumSum[d] += mass * VEL(c, i, d);
                }

                m_vSqd_Sum += mass * temp;
            }
        }

    /* update GroupMomentum & GroupParticlesNum & GroupKineticEnergy & GroupTemperature */
    updateGroupProperties(md, MomentumSum, ParticlesNum, m_vSqd_Sum);

    /* check if the particles must be placed in other cells or subdomains */
    _fmd_refreshGrid(md);
}

static void SymplecticEuler_integrate(fmd_t *md, fmd_real_t duration)
{
    fmd_real_t start = md->time;

    if (md->active_ttm_turi != NULL) _fmd_ttm_getReady(md);

    /* update force-independent fields */
    if (md->turies_num > 0)
        _fmd_turies_update(md, md->time_iteration, md->time, FMD_TRUE, FMD_TRUE, FMD_FALSE);

    /* compute forces for the first time */
    _fmd_dync_updateForces(md);

    /* update force-dependent fields */
    if (md->turies_num > 0)
        _fmd_turies_update(md, md->time_iteration, md->time, FMD_TRUE, FMD_TRUE, FMD_TRUE);

    /* make a timer-tick event */
    if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);

    while (md->time < start + duration)
    {
        /* update positions and velocities */
        SymplecticEuler_takeOneStep(md);

        /* update force-independent fields */
        if (md->turies_num > 0)
            _fmd_turies_update(md, md->time_iteration + 1, md->time + md->timestep,
                               FMD_TRUE, FMD_TRUE, FMD_FALSE);

        /* compute forces */
        _fmd_dync_updateForces(md);

        /* step completed: now update time and time_iteration */
        md->time += md->timestep;
        md->time_iteration++;

        /* update force-dependent fields */
        if (md->turies_num > 0)
            _fmd_turies_update(md, md->time_iteration, md->time, FMD_TRUE, FMD_TRUE, FMD_TRUE);

        /* make a timer-tick event */
        if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);
    }
    /* end of the time loop */
}

static void VelocityVerlet_integrate(fmd_t *md, fmd_real_t duration, fmd_bool_t UseThermostat)
{
    fmd_real_t start = md->time;

    /* compute forces for the first time */
    _fmd_dync_updateForces(md);

    /* update fields */
    if (md->turies_num > 0)
        _fmd_turies_update(md, md->time_iteration, md->time, FMD_TRUE, FMD_TRUE, FMD_TRUE);

    /* make a timer-tick event */
    if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);

    while (md->time < start + duration)
    {
        /* take first step of velocity Verlet integrator */
        VelocityVerlet_startStep(md, UseThermostat);

        /* compute forces */
        _fmd_dync_updateForces(md);

        /* take last step of velocity Verlet integrator */
        VelocityVerlet_finishStep(md);

        /* step completed: now update time and time_iteration */
        md->time += md->timestep;
        md->time_iteration++;

        /* update fields */
        if (md->turies_num > 0)
            _fmd_turies_update(md, md->time_iteration, md->time, FMD_TRUE, FMD_TRUE, FMD_TRUE);

        /* make a timer-tick event */
        if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);
    }
    /* end of the time loop */
}

static void setActiveGroup(fmd_t *md, int GroupID)
{
    if (GroupID != md->ActiveGroup && md->ActiveGroup != FMD_GROUP_ALL)
    {
        md->ActiveGroup = GroupID;
        _fmd_compute_GroupTemperature_etc_localgrid(md);
    }
    else
        md->ActiveGroup = GroupID;
}

void fmd_dync_integrate(fmd_t *md, int GroupID, fmd_real_t duration, fmd_real_t timestep)
{
    if (!md->ParticlesDistributed) _fmd_matt_distribute(md);

    setActiveGroup(md, GroupID);

    md->timestep = timestep;

    if (md->active_ttm_turi != NULL)
        SymplecticEuler_integrate(md, duration);
    else
        VelocityVerlet_integrate(md, duration, FMD_FALSE);
}

void fmd_dync_equilibrate(fmd_t *md, int GroupID, fmd_real_t duration,
  fmd_real_t timestep, fmd_real_t strength, fmd_real_t temperature)
{
    fmd_real_t bak_time;

    if (!md->ParticlesDistributed) _fmd_matt_distribute(md);

    // make backups
    bak_time = md->time;

    // initialize
    md->time = 0.0;
    md->timestep = timestep;
    md->DesiredTemperature = temperature;
    md->BerendsenThermostatParam = strength;

    setActiveGroup(md, GroupID);

    VelocityVerlet_integrate(md, duration, FMD_TRUE);

    // restore backups
    md->time = bak_time;
}
