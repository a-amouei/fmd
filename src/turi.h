/*
  turi.h: This file is part of Free Molecular Dynamics

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

#ifndef TURI_H
#define TURI_H

#include "config.h"

typedef enum
{
    FMD_FIELD_NUMBER,
    FMD_FIELD_TEMPERATURE
} fmd_field_t; /* category of the field */

typedef struct
{
    fmd_field_t cat;
    void *data;
} field_t;

typedef enum
{
    FMD_TURI_CUSTOM,
    FMD_TURI_TTM
} fmd_turi_t; /* category of the turi */

typedef struct _turi
{
    fmd_turi_t cat;
    int tdims_global[3];        /* global dimenstions of the turi */
    int tdims[3];               /* dimenstions of the turi in current subdomain.
                                   The subdomain may share some of its turi-cells with
                                   neighbor subdomains. */
    int tcell_start[3];         /* global index of the first turi-cell in current subdomain */
    int tcell_stop[3];          /* global index of the first turi-cell outside current subdomain */
    double tcellh[3];           /* size of each cell of the turi (global) */
    unsigned fields_num;
    field_t *fields;
} turi_t;

typedef struct _fmd fmd_t;

unsigned fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz);

#endif /* TURI_H */
