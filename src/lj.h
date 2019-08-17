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

#define LJ_PAIR_UPDATE_FORCE_AND_POTENERGY                                      \
    {                                                                           \
        COMPUTE_rv_AND_r2;                                                      \
                                                                                \
        LJ_6_12_t *lj = (LJ_6_12_t *)pottable[atomkind1][atomkind2].data;       \
                                                                                \
        if (r2 < lj->cutoff_sqr)                                                \
        {                                                                       \
            double inv_r2, inv_rs2, inv_rs6, inv_rs12;                          \
                                                                                \
            /* force, F = -(d/dr)U */                                           \
            inv_r2 = 1.0/r2;                                                    \
            inv_rs2 = SQR(lj->sig) * inv_r2;                                    \
            inv_rs6 = inv_rs2 * inv_rs2 * inv_rs2;                              \
            inv_rs12 = SQR(inv_rs6);                                            \
            double factor = 48.0 * lj->eps * inv_r2 * (inv_rs12 - 0.5*inv_rs6); \
            for (d=0; d<3; d++)                                                 \
                item1_p->F[d] += rv[d] * factor;                                \
                                                                                \
            /* potential energy, U = 4*eps*( (sig/r)^12 - (sig/r)^6 ) */        \
            potEnergy += 4.0 * lj->eps * (inv_rs12 - inv_rs6);                  \
        }                                                                       \
    }

typedef struct
{
    double eps;
    double sig;
    double cutoff_sqr;
} LJ_6_12_t;

typedef struct fmd_t fmd_t;

void fmd_computeLJ(fmd_t *md);

#endif /* LJ_H */
