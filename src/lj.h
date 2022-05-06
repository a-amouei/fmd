/*
  lj.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

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

#ifndef LJ_H
#define LJ_H

#include "config.h"
#include "types.h"

#define LJ_PAIR_UPDATE_FORCE_AND_POTENERGY                                          \
    {                                                                               \
        COMPUTE_rv_AND_r2;                                                          \
                                                                                    \
        LJ_6_12_t *lj = (LJ_6_12_t *)pottable[atomkind1][atomkind2].data;           \
                                                                                    \
        if (r2 < lj->cutoff_sqr)                                                    \
        {                                                                           \
            fmd_real_t inv_r2, inv_rs2, inv_rs6, inv_rs12;                          \
                                                                                    \
            /* force, F = -(d/dr)U */                                               \
            inv_r2 = 1.0/r2;                                                        \
            inv_rs2 = sqrr(lj->sig) * inv_r2;                                       \
            inv_rs6 = inv_rs2 * inv_rs2 * inv_rs2;                                  \
            inv_rs12 = sqrr(inv_rs6);                                               \
            fmd_real_t factor = 48.0 * lj->eps * inv_r2 * (inv_rs12 - 0.5*inv_rs6); \
                                                                                    \
            for (int d=0; d<3; d++)                                                 \
                FRC(c1, i1, d) += rv[d] * factor;                                   \
                                                                                    \
            /* potential energy, U = 4*eps*( (sig/r)^12 - (sig/r)^6 ) */            \
            PotEnergy += 4.0 * lj->eps * (inv_rs12 - inv_rs6);                      \
        }                                                                           \
    }

typedef struct
{
    fmd_real_t eps;
    fmd_real_t sig;
    fmd_real_t cutoff_sqr;
} LJ_6_12_t;

typedef struct _fmd fmd_t;

void fmd_computeLJ(fmd_t *md);

#endif /* LJ_H */
