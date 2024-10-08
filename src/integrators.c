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

fmd_real_t fmd_dync_getTimestep(fmd_t *md)
{
    return md->timestep;
}

fmd_real_t fmd_dync_getTime(fmd_t *md)
{
    return md->time;
}

static void VelocityVerlet_startStep(fmd_t *md, bool UseThermostat)
{
    fmd_real_t VelocityScale, mass;

    if (UseThermostat)
    {
        if (!md->KineticEnergyUpdated) _fmd_compute_GroupTemperature_etc_localgrid(md);

        VelocityScale = sqrt(1 + md->timestep / md->BerendsenThermostatParam *
                        (md->DesiredTemperature / md->GroupTemperature - 1));
    }

    cell_t *c;
    unsigned i;

    for (int ic=0; ic < md->subd.nc; ic++)
        for (c = md->subd.grid + ic, i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            for (int d=0; d<DIM; d++)
            {
                if (UseThermostat) VEL(c, i, d) *= VelocityScale;

                VEL(c, i, d) += md->timestep * 0.5 / mass * FRC(c, i, d);

                POS(c, i, d) += md->timestep * VEL(c, i, d);
            }
        }

    /* check if the particles must be placed in other cells or subdomains */
    _fmd_refreshGrid(md);
}

/* updates GroupMomentum & GroupParticlesNum & GroupKineticEnergy & GroupTemperature */
void updateGroupProperties(fmd_t *md, fmd_rtuple_t MomentumSum, int ParticlesNum, fmd_real_t m_vSqd_Sum)
{
    fmd_real_t m_vSqd_SumSum;

    MPI_Allreduce(MomentumSum, md->GroupMomentum, DIM, FMD_MPI_REAL, MPI_SUM, md->MD_comm);

    MPI_Allreduce(&ParticlesNum, &md->GroupParticlesNum, 1, MPI_INT, MPI_SUM, md->MD_comm);

    MPI_Allreduce(&m_vSqd_Sum, &m_vSqd_SumSum, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
    md->GroupKineticEnergy = 0.5 * m_vSqd_SumSum;

    md->GroupTemperature = m_vSqd_SumSum / (3.0 * md->GroupParticlesNum * K_BOLTZMANN);

    md->KineticEnergyUpdated = true;
}

static void VelocityVerlet_finishStep(fmd_t *md)
{
    fmd_real_t m_vSqd_Sum = 0;
    fmd_real_t mass;
    int ParticlesNum = 0;
    fmd_rtuple_t MomentumSum = {0., 0., 0.};
    cell_t *c;
    unsigned i;

    for (int ic=0; ic < md->subd.nc; ic++)
        for (c = md->subd.grid + ic, i=0; i < c->parts_num; i++)
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
    cell_t *c;
    unsigned i;

    for (int ic=0; ic < md->subd.nc; ic++)
        for (c = md->subd.grid + ic, i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            fmd_real_t mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            for (int d=0; d<DIM; d++)
            {
                VEL(c, i, d) += md->timestep * FRC(c, i, d) / mass;
                POS(c, i, d) += md->timestep * VEL(c, i, d);
            }
        }

    /* check if the particles must be placed in other cells or subdomains */
    _fmd_refreshGrid(md);

    /* update GroupMomentum & GroupParticlesNum & GroupKineticEnergy & GroupTemperature */

    int ParticlesNum = 0;
    fmd_rtuple_t MomentumSum = {0., 0., 0.};
    fmd_real_t m_vSqd_Sum = 0;

    for (int ic=0; ic < md->subd.nc; ic++)
        for (c = md->subd.grid + ic, i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            ParticlesNum++;

            fmd_real_t mass = md->potsys.atomkinds[c->atomkind[i]].mass;
            fmd_real_t temp = 0.0;

            for (int d=0; d<DIM; d++)
            {
                temp += sqrr(VEL(c, i, d));
                MomentumSum[d] += mass * VEL(c, i, d);
            }

            m_vSqd_Sum += mass * temp;
        }

    updateGroupProperties(md, MomentumSum, ParticlesNum, m_vSqd_Sum);
}

static void SymplecticEuler_integrate(fmd_t *md, fmd_real_t duration)
{
    fmd_real_t start = md->time;

    if (md->active_ttm_turi != NULL) _fmd_ttm_getReady(md);

    /* update force-independent fields */
    if (md->turies_num > 0)
        _fmd_turies_update(md, true, true, false);

    /* compute forces for the first time */
    _fmd_dync_updateForces(md);

    /* update force-dependent fields */
    if (md->turies_num > 0)
        _fmd_turies_update(md, false, false, true);

    /* make a timer-tick event */
    if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);

    while (md->time < start + duration)
    {
        /* update positions and velocities */
        SymplecticEuler_takeOneStep(md);

        /* update time and time_iteration */
        md->time += md->timestep;
        md->time_iteration++;

        /* update force-independent fields */
        if (md->turies_num > 0)
            _fmd_turies_update(md, true, true, false);

        /* compute forces */
        _fmd_dync_updateForces(md);

        /* update force-dependent fields */
        if (md->turies_num > 0)
            _fmd_turies_update(md, false, false, true);

        /* make a timer-tick event */
        if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);
    }
    /* end of the time loop */
}

static void VelocityVerlet_integrate(fmd_t *md, fmd_real_t duration, bool UseThermostat)
{
    fmd_real_t start = md->time;

    /* compute forces for the first time */
    _fmd_dync_updateForces(md);

    /* update fields */
    if (md->turies_num > 0)
        _fmd_turies_update(md, true, true, true);

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

        /* update time and time_iteration */
        md->time += md->timestep;
        md->time_iteration++;

        /* update fields */
        if (md->turies_num > 0)
            _fmd_turies_update(md, true, true, true);

        /* make a timer-tick event */
        if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);
    }
    /* end of the time loop */
}

static void setActiveGroup(fmd_t *md, int GroupID)
{
    if (GroupID != md->ActiveGroup)
    {
        md->ActiveGroup = GroupID;
        md->KineticEnergyUpdated = false;
    }
}

static void create_force_arrays_of_cells(fmd_t *md)
{
    bool FembP_alter = md->cellinfo.vaream_active != md->potsys.eam_applied;

    if (!md->cellinfo.F_active || FembP_alter)
        for (int ic=0; ic < md->subd.ncm; ic++)
            _fmd_cell_create_force_arrays(md, md->subd.grid + ic, FembP_alter);

    md->cellinfo.F_active = true;
    md->cellinfo.vaream_active = md->potsys.eam_applied;
}

void fmd_dync_integrate(fmd_t *md, int GroupID, fmd_real_t duration, fmd_real_t timestep)
{
    if (!md->ParticlesDistributed) _fmd_matt_distribute(md);
    _fmd_pot_update_and_process_potcats(md);
    create_force_arrays_of_cells(md);

    setActiveGroup(md, GroupID);

    md->timestep = timestep;

    if (md->active_ttm_turi != NULL)
        SymplecticEuler_integrate(md, duration);
    else
        VelocityVerlet_integrate(md, duration, false);
}

void fmd_dync_equilibrate(fmd_t *md, int GroupID, fmd_real_t duration,
  fmd_real_t timestep, fmd_real_t tau, fmd_real_t temperature)
{
    if (!md->ParticlesDistributed) _fmd_matt_distribute(md);
    _fmd_pot_update_and_process_potcats(md);
    create_force_arrays_of_cells(md);

    // make backups
    fmd_real_t bak_time = md->time;

    // initialize
    md->time = 0.0;
    md->timestep = timestep;
    md->DesiredTemperature = temperature;
    md->BerendsenThermostatParam = tau;

    setActiveGroup(md, GroupID);

    VelocityVerlet_integrate(md, duration, true);

    // restore backups
    md->time = bak_time;
}
