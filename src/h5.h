/*
  h5.h: This file is part of Free Molecular Dynamics

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

#ifndef H5_H
#define H5_H

#include "config.h"
#include <hdf5.h>
#include "types.h"

typedef struct
{
    hid_t ds_scalar;
    hid_t ds_simple_1_3;
} h5_dataspaces_t;

typedef struct _turi turi_t;
typedef struct _fmd_array3s fmd_array3s_t;
typedef struct _fmd fmd_t;

void _fmd_h5_ds_init(fmd_t *md, h5_dataspaces_t *ds);
void _fmd_h5_ds_free(fmd_t *md, h5_dataspaces_t *ds);
void _fmd_h5_save_scalar_field_float(fmd_t *md, fmd_string_t fieldname, turi_t *t, fmd_string_t path, fmd_array3s_t *arr);
void _fmd_h5_save_tuple_field_float(fmd_t *md, fmd_string_t fieldname, turi_t *t, fmd_string_t path, fmd_array3s_t *arr);

#endif /* H5_H */
