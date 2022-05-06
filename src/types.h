/*
  types.h: This file is part of Free Molecular Dynamics

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

#ifndef TYPES_H
#define TYPES_H

#include "config.h"

#define FMD_TRUE 1
#define FMD_FALSE 0

#ifdef FMD_FLOAT_CALCS
#define FMD_MPI_REAL  MPI_FLOAT
  typedef float fmd_real_t;
#else
#define FMD_MPI_REAL  MPI_DOUBLE
  typedef double fmd_real_t;
#endif

typedef char *fmd_string_t;
typedef int fmd_bool_t;
typedef void *fmd_pointer_t;

typedef int fmd_itriple_t[3];
typedef fmd_itriple_t fmd_ituple_t;

typedef unsigned fmd_utriple_t[3];
typedef fmd_utriple_t fmd_utuple_t;

typedef fmd_real_t fmd_rtriple_t[3];
typedef fmd_rtriple_t fmd_rtuple_t;

typedef fmd_bool_t fmd_btriple_t[3];
typedef fmd_btriple_t fmd_btuple_t;

typedef float fmd_ftriple_t[3];

typedef int fmd_handle_t;

typedef enum
{
    DATATYPE_REAL,
    DATATYPE_FLOAT,
    DATATYPE_DOUBLE,
    DATATYPE_INT,
    DATATYPE_UNSIGNED,
    DATATYPE_RTUPLE,
    DATATYPE_CELL             /* a cell_t */
} datatype_t;

#endif /* TYPES_H */
