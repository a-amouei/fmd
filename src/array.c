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
