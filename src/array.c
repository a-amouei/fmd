/*
  array.c: This file is part of Free Molecular Dynamics

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

#include "array.h"
#include <stdlib.h>
#include <assert.h>

/* creates an uninitialized "neat" dim1 x dim2 array;
   elsize is the size of each data element in the array in bytes;
   returns NULL if unsuccessful. */
void **_fmd_array_neat2d_create(unsigned dim1, unsigned dim2, unsigned elsize)
{
    void **arr = (void **)malloc(dim1*sizeof(void *) + dim1*dim2*elsize);

    if (arr != NULL)
    {
        // ptr points to the first data element in arr
        char *ptr = (char *)(arr + dim1);

        for (int i=0; i<dim1; i++)
            arr[i] = (void *)(ptr + i * dim2 * elsize);
    }

    return arr;
}

void _fmd_array_neat2d_free(void **array)
{
    free((void *)array);
}

fmd_pointer_t ***_fmd_array_neat3d_pointer_create(unsigned dim1, unsigned dim2, unsigned dim3)
{
    size_t s = dim1*sizeof(fmd_pointer_t **) + dim1*dim2*sizeof(fmd_pointer_t *) + dim1*dim2*dim3*sizeof(fmd_pointer_t);

    fmd_pointer_t ***arr = (fmd_pointer_t ***)malloc(s);

    if (arr != NULL)
    {
        for (int i=0; i<dim1; i++)
        {
            arr[i] = (fmd_pointer_t **)(arr + dim1 + i * (dim2 + dim2*dim3));

            for (int j=0; j<dim2; j++)
            {
                arr[i][j] = (fmd_pointer_t *)(arr[i] + dim2 + j * dim3);

                for (int k=0; k<dim3; k++)
                    arr[i][j][k] = NULL;
            }
        }
    }

    return arr;
}

void _fmd_array_neat3d_pointer_free(fmd_pointer_t ***array)
{
    free((void *)array);
}

fmd_pointer_t ***_fmd_array_ordinary3d_pointer_create(unsigned dim1, unsigned dim2, unsigned dim3)
{
    fmd_pointer_t ***arr;

    arr = (fmd_pointer_t ***)malloc(dim1 * sizeof(fmd_pointer_t **));

    if (arr != NULL)
        for (int i=0; i < dim1; i++)
        {
            arr[i] = (fmd_pointer_t **)malloc(dim2 * sizeof(fmd_pointer_t *));
            assert(arr[i] != NULL);
            /* TO-DO: handle memory error */

            for (int j=0; j < dim2; j++)
            {
                arr[i][j] = (fmd_pointer_t *)malloc(dim3 * sizeof(fmd_pointer_t));
                assert(arr[i][j] != NULL);
                /* TO-DO: handle memory error */

                for (int k=0; k < dim3; k++)
                    arr[i][j][k] = NULL;
            }
        }

    return arr;
}

void _fmd_array_ordinary3d_pointer_free(fmd_pointer_t ***array, unsigned dim1, unsigned dim2, unsigned dim3)
{
    for (int i=0; i < dim1; i++)
    {
        for (int j=0; j < dim2; j++)
            free(array[i][j]);
        free(array[i]);
    }

    free(array);
}

/* This function first tries to create a "neat" 3D pointer-array. If it isn't possible,
   then makes an ordinary one. */
fmd_pointer_t ***_fmd_array_3d_pointer_create(unsigned dim1, unsigned dim2, unsigned dim3, array_kind_t *type)
{
    fmd_pointer_t ***arr;

    arr = _fmd_array_neat3d_pointer_create(dim1, dim2, dim3);
    if (arr != NULL)
        *type = ARRAY_NEAT3D;
    else
    {
        arr = _fmd_array_ordinary3d_pointer_create(dim1, dim2, dim3);
        *type = ARRAY_ORDINARY3D;
    }

    return arr;
}

void _fmd_array_3d_pointer_free(fmd_pointer_t ***array, array_kind_t type, unsigned dim1, unsigned dim2, unsigned dim3)
{
    switch (type)
    {
        case ARRAY_NEAT3D:
            _fmd_array_neat3d_pointer_free(array);
            break;
        case ARRAY_ORDINARY3D:
            _fmd_array_ordinary3d_pointer_free(array, dim1, dim2, dim3);
            break;
    }
}
