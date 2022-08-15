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

typedef struct _ttm
{
    fmd_real_t C_gamma;     /* used when electron heat capacity is linear function of temperature */
    fmd_real_t K;           /* used when electron heat conductivity is constant */
    fmd_real_t G;           /* used when electron-ion coupling factor is constant */
} ttm_t;

typedef struct _turi turi_t;

ttm_t *_fmd_ttm_constructor(turi_t *t);
void _fmd_ttm_destructor(ttm_t **ttm);

#endif /* TTM_H */
