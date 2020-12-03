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

#define LOOP3D(iv, minv, upv)                                         \
    for ( (iv)[0]=(minv)[0]; (iv)[0]<(upv)[0]; (iv)[0]++ )            \
        for ( (iv)[1]=(minv)[1]; (iv)[1]<(upv)[1]; (iv)[1]++ )        \
            for ( (iv)[2]=(minv)[2]; (iv)[2]<(upv)[2]; (iv)[2]++ )

inline fmd_real_t sqrr(fmd_real_t x)
{
    return x*x;
}

extern const fmd_itriple_t _fmd_ThreeZeros_int;
extern const fmd_ftriple_t _fmd_ThreeZeros_float;

#endif /* GENERAL_H */
