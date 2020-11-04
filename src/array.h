/*
  array.h: This file is part of Free Molecular Dynamics

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

#ifndef ARRAY_H
#define ARRAY_H

#include "config.h"
#include "types.h"

typedef enum {ARRAY_NEAT3D, ARRAY_SEMINEAT3D, ARRAY_ORDINARY3D} array_kind_t;

typedef struct
{
    void ***data;
    unsigned dim1;
    unsigned dim2;
    unsigned dim3;
    array_kind_t kind;
} fmd_array3D_t;

void **_fmd_array_neat2d_create(unsigned dim1, unsigned dim2, unsigned elsize);
void _fmd_array_neat2d_free(void **array);
void ***_fmd_array_ordinary3d_create(unsigned dim1, unsigned dim2, unsigned dim3, unsigned elsize);
void _fmd_array_ordinary3d_free(void ***array, unsigned dim1, unsigned dim2);
void _fmd_array_3d_create(unsigned dim1, unsigned dim2, unsigned dim3, unsigned elsize, fmd_array3D_t *array);
void _fmd_array_3d_free(fmd_array3D_t *array);
void ***_fmd_array_neat3d_create(unsigned dim1, unsigned dim2, unsigned dim3, unsigned elsize);
void _fmd_array_neat3d_free(void ***array);
void _fmd_array_3d_pointer_clean(fmd_pointer_t ***array, unsigned dim1, unsigned dim2, unsigned dim3);
void _fmd_array_3d_rtuple_clean(fmd_rtuple_t ***array, unsigned dim1, unsigned dim2, unsigned dim3);
void _fmd_array_3d_real_clean(fmd_real_t ***array, unsigned dim1, unsigned dim2, unsigned dim3);
void _fmd_array_3d_unsigned_clean(unsigned ***array, unsigned dim1, unsigned dim2, unsigned dim3);
void ***_fmd_array_semineat3d_create(unsigned dim1, unsigned dim2, unsigned dim3, unsigned elsize);
void _fmd_array_semineat3d_free(void ***array, unsigned dim1);

#endif /* ARRAY_H */
