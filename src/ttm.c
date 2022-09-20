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
#include "ttm.h"
#include "general.h"
#include "turi.h"
#include "turi_ghost.h"
#include "types.h"
#include "base.h"

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

void _fmd_ttm_destructor(ttm_t **ttm)
{
    _fmd_array_3d_free(&(*ttm)->Te_aux);
    free(*ttm);
    *ttm = NULL;
}

ttm_t *_fmd_ttm_constructor(fmd_t *md, turi_t *t)
{
    ttm_t *ttm = (ttm_t *)m_alloc(sizeof(ttm_t));

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1:
        {
            int inum = _fmd_field_add(t, FMD_FIELD_NUMBER, md->timestep, FMD_FALSE);
            int ivcm = _fmd_field_add(t, FMD_FIELD_VCM, md->timestep, FMD_TRUE);
            int iTi = _fmd_field_add(t, FMD_FIELD_TEMPERATURE, md->timestep, FMD_FALSE);
            int iTe = _fmd_field_add(t, FMD_FIELD_TTM_TE, md->timestep, FMD_FALSE);
            int ixi = _fmd_field_add(t, FMD_FIELD_TTM_XI, md->timestep, FMD_TRUE);

            _fmd_array_3d_create(t->tdims_local, sizeof(fmd_real_t), DATATYPE_REAL, &ttm->Te_aux);
            assert(ttm->Te_aux.data != NULL);

            fmd_rtuple_t ll, ul;

            fmd_matt_findLimits(md, ll, ul);
            ttm->frontsurf = ll[DIM-1];

            if (t->tdims_global[0] == 1 && t->tdims_global[1] == 1)
            {
                ttm->dim = 1;

                ttm->num_1d = ((unsigned ***)t->fields[inum].data.data)[0][0];
                ttm->vcm_1d = ((fmd_rtuple_t ***)t->fields[ivcm].data.data)[0][0];
                ttm->Ti_1d = ((fmd_real_t ***)t->fields[iTi].data.data)[0][0];
                ttm->Te_1d = ((fmd_real_t ***)t->fields[iTe].data.data)[0][0];
                ttm->Te2_1d = ((fmd_real_t ***)ttm->Te_aux.data)[0][0];
                ttm->xi_1d = ((fmd_real_t ***)t->fields[ixi].data.data)[0][0];
            }
            else
            {
                ttm->dim = 3;

                ttm->num = (unsigned ***)t->fields[inum].data.data;
                ttm->vcm = (fmd_rtuple_t ***)t->fields[ivcm].data.data;
                ttm->Ti = (fmd_real_t ***)t->fields[iTi].data.data;
                ttm->Te = (fmd_real_t ***)t->fields[iTe].data.data;
                ttm->Te2 = (fmd_real_t ***)ttm->Te_aux.data;
                ttm->xi = (fmd_real_t ***)t->fields[ixi].data.data;
            }
        }
            break;

        default:
            assert(0);
    }

    return ttm;
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

    switch (cat)
    {
        case FMD_TTM_TE_CONSTANT:
        {
            fmd_ituple_t itc;

            LOOP3D(itc, t->itc_start, t->itc_stop)
                ARRAY_ELEMENT(ttm->Te, itc) = ((fmd_ttm_params_Te_constant_t *)params)->value;
        }
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

static void pack(fmd_t *md, turi_t *t, fmd_ituple_t vitc_start, fmd_ituple_t vitc_stop, size_t *size, fmd_pointer_t *out)
{
    *size = (vitc_stop[0]-vitc_start[0]) * (vitc_stop[1]-vitc_start[1]) * (vitc_stop[2]-vitc_start[2]);
    *out = malloc(*size);

    printf("process %d sends %ld turi-cell!\n", md->SubDomain.myrank, *size);
}

static void unpack(fmd_t *md, turi_t *t, fmd_ituple_t vitc_start, fmd_ituple_t vitc_stop, fmd_pointer_t in)
{

}

void ttm_dummy(fmd_t *md, fmd_handle_t ti)
{
    turi_t *t = &md->turies[ti];

    if (t->ownerscomm.owned_tcells_num > 0)
        _fmd_turi_update_ghosts(md, t, pack, unpack);
}
