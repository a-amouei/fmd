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

void _fmd_ttm_destructor(ttm_t **ttm)
{
    free(*ttm);
    *ttm = NULL;
}

ttm_t *_fmd_ttm_constructor(turi_t *t)
{
    ttm_t *ttm = (ttm_t *)m_alloc(sizeof(ttm_t));

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
