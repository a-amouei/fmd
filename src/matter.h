/*
  matter.h: This file is part of Free Molecular Dynamics

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

#ifndef MATTER_H
#define MATTER_H

#include "config.h"

#include "types.h"

typedef struct _fmd fmd_t;

typedef struct
{
    float x[DIM];
    float var;
    unsigned atomkind;
} config_atom_t;

void fmd_matt_findLimits(fmd_t *md, int GroupID, fmd_rtuple_t LowerLimit, fmd_rtuple_t UpperLimit);
void _fmd_matt_distribute(fmd_t *md);
void _fmd_compute_GroupTemperature_etc_localgrid(fmd_t *md);

#endif /* MATTER_H */
