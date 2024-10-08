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
#include "general.h"

typedef struct _fmd fmd_t;

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

void _fmd_array_3d_pointer_clean(fmd_pointer_t ***array, fmd_utriple_t dims)
{
    fmd_utriple_t iv;

    LOOP3D(iv, _fmd_ThreeZeros_int, dims)
        for (int d=0; d<3; d++)
            ARRAY_ELEMENT(array, iv) = NULL;
}

void _fmd_array_3d_rtuple_clean(fmd_rtuple_t ***array, fmd_utriple_t dims)
{
    fmd_utriple_t iv;

    LOOP3D(iv, _fmd_ThreeZeros_int, dims)
        for (int d=0; d<3; d++)
            ARRAY_ELEMENT(array, iv)[d] = 0.0;
}

void _fmd_array_3d_real_clean(fmd_real_t ***array, fmd_utriple_t dims)
{
    fmd_utriple_t iv;

    LOOP3D(iv, _fmd_ThreeZeros_int, dims)
        ARRAY_ELEMENT(array, iv) = 0.0;
}

void _fmd_array_3d_unsigned_clean(unsigned ***array, fmd_utriple_t dims)
{
    fmd_utriple_t iv;

    LOOP3D(iv, _fmd_ThreeZeros_int, dims)
        ARRAY_ELEMENT(array, iv) = 0;
}

void ***_fmd_array_neat3d_create(fmd_utriple_t dims, unsigned elsize)
{
    size_t s_ptrs1 = dims[0] * sizeof(void **);
    size_t s_ptrs2 = dims[1] * sizeof(void *);
    size_t s_2d_arr = s_ptrs2 + dims[1] * dims[2] * elsize;

    void ***arr = (void ***)malloc(s_ptrs1 + dims[0] * s_2d_arr);

    if (arr != NULL)
    {
        for (int i=0; i < dims[0]; i++)
        {
            arr[i] = (void **)((char *)arr + s_ptrs1 + i * s_2d_arr);

            for (int j=0; j < dims[1]; j++)
                arr[i][j] = (void *)((char *)arr[i] + s_ptrs2 + j * dims[2] * elsize);
        }
    }

    return arr;
}

void _fmd_array_neat3d_free(void ***array)
{
    free((void *)array);
}

void ***_fmd_array_semineat3d_create(fmd_utriple_t dims, unsigned elsize)
{
    void ***arr;

    arr = (void ***)malloc(dims[0] * sizeof(void **));

    if (arr != NULL)
    {
        for (unsigned i=0; i < dims[0]; dims[0]++)
        {
            arr[i] = _fmd_array_neat2d_create(dims[1], dims[2], elsize);

            if (arr[i] == NULL)  /* free memory and return NULL */
            {
                for (unsigned j=0; j<i; j++)
                    _fmd_array_neat2d_free(arr[i]);

                free(arr);

                return NULL;
            }
        }
    }

    return arr;
}

void _fmd_array_semineat3d_free(void ***array, unsigned dim1)
{
    for (unsigned i=0; i < dim1; i++)
        _fmd_array_neat2d_free(array[i]);
    free((void *)array);
}

void ***_fmd_array_ordinary3d_create(fmd_t *md, fmd_utriple_t dims, unsigned elsize)
{
    void ***arr;

    arr = (void ***)m_alloc(md, dims[0] * sizeof(void **));

    for (int i=0; i < dims[0]; i++)
    {
        arr[i] = (void **)m_alloc(md, dims[1] * sizeof(void *));

        for (int j=0; j < dims[1]; j++)
            arr[i][j] = (void *)m_alloc(md, dims[2] * elsize);
    }

    return arr;
}

void _fmd_array_ordinary3d_free(void ***array, unsigned dim1, unsigned dim2)
{
    for (int i=0; i < dim1; i++)
    {
        for (int j=0; j < dim2; j++)
            free(array[i][j]);
        free(array[i]);
    }

    free(array);
}

/* This function first tries to create a "neat" 3D array. If it isn't possible,
   then tries to make a "semi-neat" one. If not possible yet, makes an "ordinary" 3D array. */
void _fmd_array_3d_create(fmd_t *md, fmd_utriple_t dims, unsigned elsize, datatype_t dt, fmd_array3s_t *array)
{
    array->data = _fmd_array_neat3d_create(dims, elsize);
    if (array->data != NULL)
        array->kind = ARRAY_NEAT3D;
    else
    {
        array->data = _fmd_array_semineat3d_create(dims, elsize);
        if (array->data != NULL)
            array->kind = ARRAY_SEMINEAT3D;
        else
        {
            array->data = _fmd_array_ordinary3d_create(md, dims, elsize);
            array->kind = ARRAY_ORDINARY3D;
        }
    }

    array->dims[0] = dims[0];
    array->dims[1] = dims[1];
    array->dims[2] = dims[2];

    array->datatype = dt;
}

void _fmd_array_3d_free(fmd_array3s_t *array)
{
    switch (array->kind)
    {
        case ARRAY_NEAT3D:
            _fmd_array_neat3d_free(array->data);
            break;
        case ARRAY_SEMINEAT3D:
            _fmd_array_semineat3d_free(array->data, array->dims[0]);
            break;
        case ARRAY_ORDINARY3D:
            _fmd_array_ordinary3d_free(array->data, array->dims[0], array->dims[1]);
            break;
    }
}

void fmd_array3s_free(fmd_array3s_t *array)
{
    if (array != NULL)
    {
        _fmd_array_3d_free(array);
        free(array);
    }
}

float *_fmd_array_convert_numerical_scalar_3d_to_flat_float(fmd_t *md, fmd_array3s_t *array)
{
    float *f = malloc(array->dims[0] * array->dims[1] * array->dims[2] * sizeof(*f));
    if (f == NULL) return NULL;

    int cc = 0;
    fmd_utriple_t iv;

    switch(array->datatype)
    {
        case DATATYPE_REAL:
            LOOP3D(iv, _fmd_ThreeZeros_int, array->dims)
                f[cc++] = ARRAY_ELEMENT((fmd_real_t ***)(array->data), iv);
            break;

        case DATATYPE_FLOAT:
            LOOP3D(iv, _fmd_ThreeZeros_int, array->dims)
                f[cc++] = ARRAY_ELEMENT((float ***)(array->data), iv);
            break;

        case DATATYPE_DOUBLE:
            LOOP3D(iv, _fmd_ThreeZeros_int, array->dims)
                f[cc++] = ARRAY_ELEMENT((double ***)(array->data), iv);
            break;

        case DATATYPE_INT:
            LOOP3D(iv, _fmd_ThreeZeros_int, array->dims)
                f[cc++] = ARRAY_ELEMENT((int ***)(array->data), iv);
            break;

        case DATATYPE_UNSIGNED:
            LOOP3D(iv, _fmd_ThreeZeros_int, array->dims)
                f[cc++] = ARRAY_ELEMENT((unsigned ***)(array->data), iv);
            break;

        default:
            _fmd_error_unacceptable_int_value(md, true, __FILE__, (fmd_string_t)__func__, __LINE__,
                                              "datatype of array elements", array->datatype);
            return NULL;
    }

    return f;
}

float *_fmd_array_convert_numerical_tuple_3d_to_flat_float(fmd_t *md, fmd_array3s_t *array)
{
    float *f = (float *)malloc(3 * array->dims[0] * array->dims[1] * array->dims[2] * sizeof(*f));
    if (f == NULL) return NULL;

    int cc = 0;
    fmd_utriple_t iv;

    switch(array->datatype)
    {
        case DATATYPE_RTUPLE:
            LOOP3D(iv, _fmd_ThreeZeros_int, array->dims)
                for (int d=0; d<DIM; d++)
                    f[cc++] = ARRAY_ELEMENT((fmd_rtuple_t ***)(array->data), iv)[d];
            break;

        default:
            _fmd_error_unacceptable_int_value(md, true, __FILE__, (fmd_string_t)__func__, __LINE__,
                                              "datatype of array elements", array->datatype);
            return NULL;
    }

    return f;
}
