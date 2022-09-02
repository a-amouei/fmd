/*
  forces.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani

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

#ifndef FORCES_H
#define FORCES_H

#include "config.h"

#define COMPUTE_rv_AND_r2                                                    \
    do                                                                       \
    {                                                                        \
        for (int d=0; d<DIM; d++)                                            \
        {                                                                    \
            if (md->ns[d] == 1)                                              \
            {                                                                \
                if (kc[d]==-1)                                               \
                    rv[d] = x1[d] - x2[d] + md->l[d];                        \
                else                                                         \
                    if (kc[d] == md->nc[d])                                  \
                        rv[d] = x1[d] - x2[d] - md->l[d];                    \
                    else                                                     \
                        rv[d] = x1[d] - x2[d];                               \
            }                                                                \
            else                                                             \
                rv[d] = x1[d] - x2[d];                                       \
        }                                                                    \
        r2 = sqrr(rv[0])+sqrr(rv[1])+sqrr(rv[2]);                            \
    } while (0)

typedef struct _fmd fmd_t;

void _fmd_dync_updateForces(fmd_t *md);

#endif /* FORCES_H */
