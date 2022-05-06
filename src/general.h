/*
  general.h: This file is part of Free Molecular Dynamics

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

#ifndef GENERAL_H
#define GENERAL_H

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "config.h"
#include "types.h"

#define RANK0 0

/* error codes */
#define ERROR_NCELL_TOO_SMALL                   1
#define ERROR_UNEXPECTED_PARTICLE_POSITION      2
#define ERROR_UNABLE_OPEN_FILE                  3
#define ERROR_UNSUITABLE_FILE                   4

#define LOOP3D(i3, min3, up3)                                         \
    for ( (i3)[0]=(min3)[0]; (i3)[0]<(up3)[0]; (i3)[0]++ )            \
        for ( (i3)[1]=(min3)[1]; (i3)[1]<(up3)[1]; (i3)[1]++ )        \
            for ( (i3)[2]=(min3)[2]; (i3)[2]<(up3)[2]; (i3)[2]++ )

/* converts index from 3D space to flat space */
#define INDEX_FLAT(i3, n3)   ((i3)[2] + (n3)[2]*((i3)[1] + (n3)[1]*(i3)[0]))

/* converts index from flat space to 3D space */
#define INDEX_3D(i, n3, i3)                                            \
    ((i3)[2]= (i)%(n3)[2],                                             \
     (i3)[1]=((i)/(n3)[2])%(n3)[1],                                    \
     (i3)[0]=((i)/(n3)[2])/(n3)[1])

extern const fmd_itriple_t _fmd_ThreeZeros_int;
extern const fmd_ftriple_t _fmd_ThreeZeros_float;

FILE *f_open(char *filename, char *modes);

static inline fmd_real_t sqrr(fmd_real_t x)
{
    return x*x;
}

static inline void *re_alloc(void *ptr, size_t size)
{
    void *res;

    res = realloc(ptr, size);
    assert(res != NULL);   /* TO-DO: handle memory error */

    return res;
}

static inline void *m_alloc(size_t size)
{
    void *res;

    res = malloc(size);
    assert(res != NULL);   /* TO-DO: handle memory error */

    return res;
}

static inline void *c_alloc(size_t nmemb, size_t size)
{
    void *res;

    res = calloc(nmemb, size);
    assert(res != NULL);   /* TO-DO: handle memory error */

    return res;
}

#endif /* GENERAL_H */
