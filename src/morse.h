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

#define MORSE_PAIR_UPDATE_FORCE_AND_POTENERGY                                               \
    {                                                                                       \
        COMPUTE_rv_AND_r2;                                                                  \
                                                                                            \
        morse_t *morse = (morse_t *)pottable[atomkind1][atomkind2].data;                    \
                                                                                            \
        if (r2 < morse->cutoff_sqr)                                                         \
        {                                                                                   \
            /* force, F = -(d/dr)U */                                                       \
            double r = sqrt(r2);                                                            \
            double inv_r = 1.0/r;                                                           \
            double exp1 = exp( -morse->alpha * (r - morse->r0) );                           \
            double exp2 = SQR(exp1);                                                        \
            double factor = 2.0 * morse->alpha * morse->D0 * inv_r * (exp2 - exp1);         \
            for (d=0; d<3; d++)                                                             \
                item1_p->F[d] += factor * rv[d];                                            \
                                                                                            \
            /* potential energy, U = D0 * ( exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)) ) */\
            potEnergy += morse->D0 * (exp2 - 2.0 * exp1);                                   \
        }                                                                                   \
    }

typedef struct
{
    double D0;
    double alpha;
    double r0;
    double cutoff_sqr;
} morse_t;

typedef struct fmd_t fmd_t;

void fmd_computeMorse(fmd_t *md);

#endif /* MORSE_H */
