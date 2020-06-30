/*
  turi.c: This file is part of Free Molecular Dynamics

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

#include "turi.h"
#include "base.h"

unsigned fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz)
{
    int i = md->turies_num;

    md->turies = (turi_t *)realloc(md->turies, (i+1) * sizeof(turi_t));
    // TO-DO: handle memory error
    assert(md->turies != NULL);

    md->turies[i].cat = cat;
    switch (cat)
    {
        case FMD_TURI_CUSTOM:
            md->turies[i].fields = NULL;
            md->turies[i].fields_num = 0;
            break;
    }

    md->turies[i].tdims_global[0] = dimx;
    md->turies[i].tdims_global[1] = dimy;
    md->turies[i].tdims_global[2] = dimz;

    md->turies_num++;

    return i;
}