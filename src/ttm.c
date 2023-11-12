/*
  ttm.c: This file is part of Free Molecular Dynamics

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

#include <stdio.h>
#include <tgmath.h>
#include "ttm.h"
#include "fmd-private.h"
#include "matter.h"
#include "general.h"
#include "turi.h"
#include "turi-ghost.h"
#include "types.h"
//#include "misc.h"

typedef struct
{
    fmd_real_t gamma;
} fmd_ttm_params_heat_capacity_linear_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_params_heat_conductivity_constant_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_params_coupling_factor_constant_t;

typedef struct
{
    unsigned value;
} fmd_ttm_params_timestep_ratio_constant_t;

typedef struct
{
    fmd_real_t fluence;
    fmd_real_t reflectance;
    fmd_real_t t0;
    fmd_real_t duration;
    fmd_real_t AbsorptionDepth;
} fmd_ttm_params_laser_simple_t;

typedef enum
{
    FMD_TTM_TE_CONSTANT
} fmd_ttm_Te_t;

typedef enum
{
    FMD_TTM_TIMESTEP_RATIO_CONSTANT
} fmd_ttm_timestep_ratio_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_params_Te_constant_t;

void _fmd_ttm_getReady(fmd_t *md)
{
    turi_t *t = md->active_ttm_turi;
    ttm_t *ttm = t->ttm;

    fmd_rtuple_t ll, ul;

    fmd_matt_findLimits(md, ll, ul);
    ttm->frontsurf = ll[DIM-1];

    fmd_real_t vol = 1.0;

    for (int d=0; d<DIM; d++)
        vol *= ul[d] - ll[d];

    ttm->initial_atoms_num = md->GroupParticlesNum / vol * t->tcell_volume;
    ttm->min_atoms_num = ttm->initial_atoms_num * ttm->CellActivFrac;
}

static void type1_1d_pack(fmd_t *md, turi_t *t, bool SendToUp, fmd_pointer_t sendbuf, int *pos)
{
    ttm_t *ttm = t->ttm;
    *pos = 0;

    int i = SendToUp ? t->itc_stop[2] - 1 : t->itc_start_owned[2];

    MPI_Pack(&ttm->num_1d[i], 1, MPI_UNSIGNED, sendbuf, INT_MAX, pos, md->MD_comm);
    MPI_Pack(&ttm->Te_1d[i],  1, FMD_MPI_REAL, sendbuf, INT_MAX, pos, md->MD_comm);
}

static void type1_1d_unpack(fmd_t *md, turi_t *t, bool SendToUp, fmd_pointer_t recvbuf)
{
    ttm_t *ttm = t->ttm;
    int pos = 0;

    int i = SendToUp ? 0 : t->itc_stop[2];

    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->num_1d[i], 1, MPI_UNSIGNED, md->MD_comm);
    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->Te_1d[i],  1, FMD_MPI_REAL, md->MD_comm);
}

static size_t type1_1d_calc_tghost_buffsize(fmd_t *md)
{
    size_t size;
    int i;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &i);
    size = i;
    MPI_Pack_size(1, FMD_MPI_REAL, md->MD_comm, &i);
    size += i;

    return size;
}

static void ttm_presolve_1d(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    if (t->has_upper_lower_owner_procs[2])
        _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp);

    /* xi intialization, cell-activation, cell-deactivation */
    for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
    {
        ttm->xi_1d[i] = 0.0;

        /* cell activation */
        if (ttm->Te_1d[i] < 0.0 && ttm->num_1d[i] >= ttm->min_atoms_num)
        {
            int count = 0;
            fmd_real_t sum = 0.0;

            int im1 = (i == t->itc_start_owned[2]) ? 0 : i-1;

            if (ttm->Te_1d[im1] >= 0.0)
            {
                sum += ttm->Te_1d[im1] * ttm->num_1d[im1];
                count += ttm->num_1d[im1];
            }

            if (ttm->Te_1d[i+1] >= 0.0)
            {
                sum += ttm->Te_1d[i+1] * ttm->num_1d[i+1];
                count += ttm->num_1d[i+1];
            }

            ttm->Te_1d[i] = (count > 0) ? sum/count : ttm->Ti_1d[i];
        }

        /* cell deactivation */
        if (ttm->num_1d[i] < ttm->min_atoms_num) ttm->Te_1d[i] = ttm->Te2_1d[i] = -1.0;
    }
}

static void ttm_type1_solve_1d(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    bool CalcSource = fabs(md->time - ttm->laser_t0) < ttm->laser_tdiff ? true : false;

    /* time loop */
    for (int j=0; j < ttm->timestep_ratio; j++)
    {
        fmd_real_t source_spcind;

        if (CalcSource)
            source_spcind = ttm->laser_factor_constant *
                            exp(ttm->laser_m_2sig2_inv * sqrr(md->time + j * ttm->timestep - ttm->laser_t0));

        /* spatial loop */
        for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
        {
            if (ttm->Te_1d[i] < 0) continue;  /* nothing to be done with an deactivated cell */

            fmd_real_t source;

            if (CalcSource)
            {
                fmd_real_t z = (i - t->itc_glob_to_loc[2] + 0.5) * t->tcellh[2];
                source = source_spcind * exp(ttm->laser_m_absdepth_inv * (z - ttm->frontsurf));
            }
            else
                source = 0.0;

            int im1 = (i == t->itc_start_owned[2]) ? 0 : i-1;
            int ilo = (ttm->num_1d[im1] < ttm->min_atoms_num) ? i : im1;
            int ihi = (ttm->num_1d[i+1] < ttm->min_atoms_num) ? i : i+1;

            ttm->Te2_1d[i] = ttm->Te_1d[i] + ttm->timestep/(ttm->C_gamma * ttm->Te_1d[i]) * (
              ttm->K * (ttm->Te_1d[ihi]-2*ttm->Te_1d[i]+ttm->Te_1d[ilo])/ttm->dz2
              - ttm->G * (ttm->Te_1d[i] - ttm->Ti_1d[i]) + source );

            ttm->xi_1d[i] += ttm->Te2_1d[i];
        } /* end of spatial loop */

        fmd_real_t *tempo = ttm->Te2_1d;
        ttm->Te2_1d = ttm->Te_1d;
        ttm->Te_1d = tempo;

        if (t->has_upper_lower_owner_procs[2] && j < ttm->timestep_ratio-1)
            _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp);
    } /* end of time loop */

    if (ttm->Te2_1d != ((fmd_real_t ***)ttm->Te_aux.data)[0][0])
    {
        fmd_array3s_t tempo = ttm->Te_aux;
        ttm->Te_aux = t->fields[ttm->iTe].data;
        t->fields[ttm->iTe].data = tempo;
    }

    for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
    {
        if (ttm->num_1d[i] >= ttm->min_atoms_num)
            ttm->xi_1d[i] = ttm->G * t->tcell_volume *
                            (ttm->xi_1d[i] / ttm->timestep_ratio - ttm->Ti_1d[i]) /
                            (3 * K_BOLTZMANN * ttm->num_1d[i] * ttm->Ti_1d[i]);
    }
}

static void ttm_init_type1(fmd_t *md, ttm_t *ttm, turi_t *t)
{
    int inum = _fmd_field_add(t, FMD_FIELD_NUMBER, md->timestep, false);
    int ivcm = _fmd_field_add(t, FMD_FIELD_VCM, md->timestep, true);
    int iTi = _fmd_field_add(t, FMD_FIELD_TEMPERATURE, md->timestep, false);
    ttm->iTe = _fmd_field_add(t, FMD_FIELD_TTM_TE, md->timestep, false);
    ttm->ixi = _fmd_field_add(t, FMD_FIELD_TTM_XI, md->timestep, true);

    _fmd_array_3d_create(t->tdims_local, sizeof(fmd_real_t), DATATYPE_REAL, &ttm->Te_aux);
    assert(ttm->Te_aux.data != NULL);

    /* set default value for "cell activation fraction" */
    ttm->CellActivFrac = 0.1;

    if (t->tdims_global[0] == 1 && t->tdims_global[1] == 1)
    {
        ttm->dim = 1;

        ttm->num_1d = ((unsigned ***)t->fields[inum].data.data)[0][0];
        ttm->vcm_1d = ((fmd_rtuple_t ***)t->fields[ivcm].data.data)[0][0];
        ttm->Ti_1d = ((fmd_real_t ***)t->fields[iTi].data.data)[0][0];
        ttm->Te_1d = ((fmd_real_t ***)t->fields[ttm->iTe].data.data)[0][0];
        ttm->Te2_1d = ((fmd_real_t ***)ttm->Te_aux.data)[0][0];
        ttm->xi_1d = ((fmd_real_t ***)t->fields[ttm->ixi].data.data)[0][0];

        ttm->num_1d[0] = ttm->num_1d[t->itc_stop[2]] = 0;

        ttm->dz2 = t->tcellh[2] * t->tcellh[2];

        ttm->preupdate_xe_te = ttm_presolve_1d;
        ttm->update_xe_te = ttm_type1_solve_1d;

       if (t->ownerscomm.owned_tcells_num == 0)
            ttm->tgp = NULL;
        else
        {
            ttm->tgp = (tghost_pack_t *)m_alloc(sizeof(tghost_pack_t));
            ttm->tgp->pack._1D = type1_1d_pack;
            ttm->tgp->unpack._1D = type1_1d_unpack;
            ttm->tgp->bufsize = type1_1d_calc_tghost_buffsize(md);
            ttm->tgp->sendbuf = m_alloc(ttm->tgp->bufsize);
            ttm->tgp->recvbuf = m_alloc(ttm->tgp->bufsize);
        }
    }
    else
    {
        ttm->dim = 3;

        ttm->num = (unsigned ***)t->fields[inum].data.data;
        ttm->vcm = (fmd_rtuple_t ***)t->fields[ivcm].data.data;
        ttm->Ti = (fmd_real_t ***)t->fields[iTi].data.data;
        ttm->Te = (fmd_real_t ***)t->fields[ttm->iTe].data.data;
        ttm->Te2 = (fmd_real_t ***)ttm->Te_aux.data;
        ttm->xi = (fmd_real_t ***)t->fields[ttm->ixi].data.data;

        ttm->preupdate_xe_te = NULL;
        ttm->update_xe_te = NULL; /* no 3D preupdater and updater is written */
    }
}

ttm_t *_fmd_ttm_construct(fmd_t *md, turi_t *t)
{
    ttm_t *ttm = (ttm_t *)m_alloc(sizeof(ttm_t));

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1:
            ttm_init_type1(md, ttm, t);
            break;

        default:
            assert(0);
    }

    return ttm;
}

void _fmd_ttm_destruct(ttm_t **ttm)
{
    _fmd_array_3d_free(&(*ttm)->Te_aux);

    if ((*ttm)->tgp != NULL)
    {
        free((*ttm)->tgp->recvbuf);
        free((*ttm)->tgp->sendbuf);
        free((*ttm)->tgp);
    }

    free(*ttm);
    *ttm = NULL;
}

void fmd_ttm_setHeatCapacity(fmd_t *md, fmd_handle_t turi, fmd_params_t *params)
{
    turi_t *t = &md->turies[turi];

    ttm_t *ttm = t->ttm;

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1:
            ttm->C_gamma = ((fmd_ttm_params_heat_capacity_linear_t *)params)->gamma * JOULE_PER_METER3_KELVIN2;
            break;

        default:
            assert(0); /* TO-DO: handle error */
    }
}

void fmd_ttm_setHeatConductivity(fmd_t *md, fmd_handle_t turi, fmd_params_t *params)
{
    turi_t *t = &md->turies[turi];

    ttm_t *ttm = t->ttm;

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1:
            ttm->K = ((fmd_ttm_params_heat_conductivity_constant_t *)params)->value * WATT_PER_METER_KELVIN;
            break;

        default:
            assert(0); /* TO-DO: handle error */
    }
}

void fmd_ttm_setCouplingFactor(fmd_t *md, fmd_handle_t turi, fmd_params_t *params)
{
    turi_t *t = &md->turies[turi];

    ttm_t *ttm = t->ttm;

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1:
            ttm->G = ((fmd_ttm_params_coupling_factor_constant_t *)params)->value * WATT_PER_METER3_KELVIN;
            break;

        default:
            assert(0); /* TO-DO: handle error */
    }
}

void fmd_ttm_setElectronTemperature(fmd_t *md, fmd_handle_t turi, fmd_ttm_Te_t cat, fmd_params_t *params)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    fmd_ituple_t itc;

    switch (cat)
    {
        case FMD_TTM_TE_CONSTANT:
            LOOP3D(itc, t->itc_start, t->itc_stop)
                ARRAY_ELEMENT((fmd_real_t ***)t->fields[ttm->iTe].data.data, itc) =
                  ((fmd_ttm_params_Te_constant_t *)params)->value;

            break;

        default:
            assert(0); /* TO-DO: handle error */
    }
}

void fmd_ttm_setTimestepRatio(fmd_t *md, fmd_handle_t turi, fmd_ttm_timestep_ratio_t cat, fmd_params_t *params)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    switch (cat)
    {
        case FMD_TTM_TIMESTEP_RATIO_CONSTANT:
            ttm->timestep_ratio = ((fmd_ttm_params_timestep_ratio_constant_t *)params)->value;
            ttm->timestep = md->timestep / ttm->timestep_ratio;
            break;

        default:
            assert(0); /* TO-DO: handle error */
    }
}

void fmd_ttm_setCellActivationFraction(fmd_t *md, fmd_handle_t turi, fmd_real_t value)
{
    assert(value >= 0 && value <= 1); /* TO-DO: handle error */

    ttm_t *ttm = md->turies[turi].ttm;

    ttm->CellActivFrac = value;
}

void fmd_ttm_setLaserSource(fmd_t *md, fmd_handle_t turi, fmd_params_t *params)
{
    turi_t *t = &md->turies[turi];

    ttm_t *ttm = t->ttm;

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1: ;
            fmd_real_t Lp = ((fmd_ttm_params_laser_simple_t *)params)->AbsorptionDepth * METER;

            fmd_real_t R = ((fmd_ttm_params_laser_simple_t *)params)->reflectance;

            fmd_real_t I0 = sqrt(4*log(2)/M_PI) * ((fmd_ttm_params_laser_simple_t *)params)->fluence /
                            ((fmd_ttm_params_laser_simple_t *)params)->duration * WATT_PER_METER2;

            fmd_real_t sigma = ((fmd_ttm_params_laser_simple_t *)params)->duration / sqrt(8*log(2)) * SECOND;

            ttm->laser_t0 = ((fmd_ttm_params_laser_simple_t *)params)->t0 * SECOND;
            ttm->laser_tdiff = 4*sqrt(2) * sigma;
            ttm->laser_m_2sig2_inv = -0.5 / sqrr(sigma);
            ttm->laser_m_absdepth_inv = -1. / Lp;
            ttm->laser_factor_constant = I0 * (1 - R) / Lp;
            break;

        default:
            assert(0); /* TO-DO: handle error */
    }
}
