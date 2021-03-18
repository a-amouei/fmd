/*
  general.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2020 Arham Amouye Foumani

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

#ifndef GENERAL_H
#define GENERAL_H

#include "config.h"
#include "types.h"

#define RANK0 0

#define LOOP3D(i3, min3, up3)                                         \
    for ( (i3)[0]=(min3)[0]; (i3)[0]<(up3)[0]; (i3)[0]++ )            \
        for ( (i3)[1]=(min3)[1]; (i3)[1]<(up3)[1]; (i3)[1]++ )        \
            for ( (i3)[2]=(min3)[2]; (i3)[2]<(up3)[2]; (i3)[2]++ )

/* converts index from 3D space to flat space */
#define INDEX_FLAT(i3, n3)   ((i3)[2] + (n3)[2]*((i3)[1] + (n3)[1]*(i3)[0]))

/* converts index from flat space to 3D space */
#define INDEX_3D(i, n3, i3)                                            \
    ((i3)[2]= (i)%(n3)[2],                                             \
     (i3)[1]=((i)/(n3)[2])%(n3)[1],                                    \
     (i3)[0]=((i)/(n3)[2])/(n3)[1])

inline fmd_real_t sqrr(fmd_real_t x)
{
    return x*x;
}

extern const fmd_itriple_t _fmd_ThreeZeros_int;
extern const fmd_ftriple_t _fmd_ThreeZeros_float;

#endif /* GENERAL_H */
