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
    for (d=0; d<3; d++)                                                      \
    {                                                                        \
        if (md->ns[d] == 1)                                                  \
        {                                                                    \
            if (kc[d]==-1)                                                   \
                rv[d] = item1_p->P.x[d] - item2_p->P.x[d] + md->l[d];        \
            else                                                             \
                if (kc[d] == md->nc[d])                                      \
                    rv[d] = item1_p->P.x[d] - item2_p->P.x[d] - md->l[d];    \
                else                                                         \
                    rv[d] = item1_p->P.x[d] - item2_p->P.x[d];               \
        }                                                                    \
        else                                                                 \
            rv[d] = item1_p->P.x[d] - item2_p->P.x[d];                       \
    }                                                                        \
    r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);

typedef struct _fmd fmd_t;

void fmd_dync_updateForces(fmd_t *md);

#endif /* FORCES_H */
