/*
  morse.h: This file is part of Free Molecular Dynamics

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

#ifndef MORSE_H
#define MORSE_H

#include "config.h"
#include "types.h"

#define MORSE_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, c1, i1, atomkind2, c2,         \
                                              i2, PotEn, pottable)                          \
    do                                                                                      \
    {                                                                                       \
        fmd_real_t r2;                                                                      \
        fmd_rtuple_t rv;                                                                    \
                                                                                            \
        fmd_real_t *x2 = &POS(c2, i2, 0);                                                   \
                                                                                            \
        COMPUTE_rv_AND_r2(x1, x2, rv, r2);                                                  \
                                                                                            \
        morse_t *morse = (morse_t *)pottable[atomkind1][atomkind2].data;                    \
                                                                                            \
        if (r2 < morse->cutoff_sqr)                                                         \
        {                                                                                   \
            /* force, F = -(d/dr)U */                                                       \
            fmd_real_t r = sqrt(r2);                                                        \
            fmd_real_t inv_r = 1.0/r;                                                       \
            fmd_real_t exp1 = exp( -morse->alpha * (r - morse->r0) );                       \
            fmd_real_t exp2 = sqrr(exp1);                                                   \
            fmd_real_t factor = 2.0 * morse->alpha * morse->D0 * inv_r * (exp2 - exp1);     \
                                                                                            \
            for (int d=0; d<DIM; d++)                                                       \
            {                                                                               \
                fmd_real_t tmp = factor * rv[d];                                            \
                FRC(c1, i1, d) += tmp;                                                      \
                FRC(c2, i2, d) -= tmp;                                                      \
            }                                                                               \
                                                                                            \
            /* potential energy, U = D0 * ( exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)) ) */\
            PotEn += morse->D0 * (exp2 - 2.0 * exp1);                                       \
        }                                                                                   \
    } while (0)

typedef struct
{
    fmd_real_t D0;
    fmd_real_t alpha;
    fmd_real_t r0;
    fmd_real_t cutoff_sqr;
} morse_t;

typedef struct _fmd fmd_t;

void _fmd_computeMorse(fmd_t *md);

#endif /* MORSE_H */
