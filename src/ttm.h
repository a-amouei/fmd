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

typedef void (*xi_te_preupdater_t)(fmd_t *md, turi_t *t, ttm_t *ttm);
typedef void (*xi_te_updater_t)(fmd_t *md, turi_t *t, ttm_t *ttm);

typedef struct
{
    fmd_real_t v2;
    fmd_real_t A;
    fmd_real_t B;
} heat_conductivity_zhigilei_t;

typedef struct _ttm
{
    int dim;                    /* 1D ttm or 3D ttm? */
    fmd_real_t C_gamma;         /* used when electron heat capacity is linear function of temperature */
    fmd_real_t *T_C;            /* used when electron heat capacity is tabulated */
    int nC;                     /* used when electron heat capacity is tabulated */
    fmd_real_t *Ct;             /* used when electron heat capacity is tabulated */
    fmd_real_t *C_DD;           /* used when electron heat capacity is tabulated */
    fmd_real_t K;               /* used when electron heat conductivity is constant */
    heat_conductivity_zhigilei_t Kzh;
    fmd_real_t G;               /* used when electron-ion coupling factor is constant */
    fmd_real_t *T_G;            /* used when electron-ion coupling factor is tabulated */
    int nG;                     /* used when electron-ion coupling factor is tabulated */
    fmd_real_t *Gt;             /* used when electron-ion coupling factor is tabulated */
    fmd_real_t *G_DD;           /* used when electron-ion coupling factor is tabulated */
    fmd_array3s_t Te_aux;
    fmd_array3s_t Kel;          /* electron heat conductivity, used when Ke changes with position */
    fmd_array3s_t Cel;          /* electron heat capacity, used when Ce changes with position */
    fmd_array3s_t Geis;         /* electron-ion coupling factor, used when G changes with position */
    int iTe;                    /* index of the field for electron temperature */
    int ixi;                    /* index of the field for xi */
    unsigned ***num;
    fmd_rtuple_t ***vcm;
    fmd_real_t ***Ti;
    fmd_real_t ***Te;
    fmd_real_t ***Te2;
    fmd_real_t ***xi;           /* electron-ion coupling factor */
    fmd_real_t ***Ke;
    fmd_real_t ***Ce;
    fmd_real_t ***Gei;
    unsigned *num_1d;
    fmd_rtuple_t *vcm_1d;
    fmd_real_t *Ti_1d;
    fmd_real_t *Te_1d;
    fmd_real_t *Te2_1d;
    fmd_real_t *xi_1d;
    fmd_real_t *Ke_1d;
    fmd_real_t *Ce_1d;
    fmd_real_t *G_1d;
    unsigned timestep_ratio;    /* the ratio of MD timestep to TTM timestep */
    fmd_real_t frontsurf;       /* position of front surface */
    fmd_real_t CellActivFrac;   /* "cell activation fraction" */
    unsigned initial_atoms_num;
    unsigned min_atoms_num;     /* ttm-cells with fewer atoms than this are deactivated */
    fmd_real_t dz2;             /* delta_z^2 -- for 1D case */
    xi_te_preupdater_t preupdate_xe_te;
    xi_te_updater_t update_xe_te;
    tghost_pack_t *tgp;
    fmd_real_t laser_factor_constant;
    fmd_real_t laser_m_absdepth_inv;  /* -1 x inverse of absorption depth */
    fmd_real_t laser_m_2sig2_inv;
    fmd_real_t laser_t0;
    fmd_real_t laser_tdiff;
} ttm_t;

ttm_t *_fmd_ttm_construct(fmd_t *md, turi_t *t);
void _fmd_ttm_destruct(turi_t *t);
void _fmd_ttm_getReady(fmd_t *md);

#endif /* TTM_H */
