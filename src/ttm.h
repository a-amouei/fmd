/*
  ttm.h: This file is part of Free Molecular Dynamics

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

#ifndef TTM_H
#define TTM_H

#include "config.h"
#include "types.h"
#include "array.h"

typedef struct _turi turi_t;
typedef struct _fmd fmd_t;
typedef struct _ttm ttm_t;
typedef struct _tghost_pack tghost_pack_t;

typedef void (*xi_te_updater_t)(fmd_t *md, turi_t *t, ttm_t *ttm);

typedef struct _ttm
{
    int dim;                    /* 1D ttm or 3D ttm? */
    fmd_real_t C_gamma;         /* used when electron heat capacity is linear function of temperature */
    fmd_real_t K;               /* used when electron heat conductivity is constant */
    fmd_real_t G;               /* used when electron-ion coupling factor is constant */
    fmd_array3D_t Te_aux;
    int iTe;                    /* index of the field for electron temperature */
    int ixi;                    /* index of the field for xi */
    unsigned ***num;
    fmd_rtuple_t ***vcm;
    fmd_real_t ***Ti;
    fmd_real_t ***Te;
    fmd_real_t ***Te2;
    fmd_real_t ***xi;           /* electron-ion coupling factor */
    unsigned *num_1d;
    fmd_rtuple_t *vcm_1d;
    fmd_real_t *Ti_1d;
    fmd_real_t *Te_1d;
    fmd_real_t *Te2_1d;
    fmd_real_t *xi_1d;
    unsigned timestep_ratio;    /* the ratio of MD timestep to TTM timestep */
    fmd_real_t timestep;
    fmd_real_t frontsurf;       /* position of front surface */
    fmd_real_t CellActivFrac;   /* "cell activation fraction" */
    unsigned initial_atoms_num;
    unsigned min_atoms_num;     /* ttm-cells with fewer atoms than this are deactivated */
    fmd_real_t dz2;             /* delta_z^2 -- for 1D case */
    xi_te_updater_t update_xe_te;
    tghost_pack_t *tgp;
} ttm_t;

ttm_t *_fmd_ttm_construct(fmd_t *md, turi_t *t);
void _fmd_ttm_destruct(ttm_t **ttm);
void _fmd_ttm_getReady(fmd_t *md);

#endif /* TTM_H */
