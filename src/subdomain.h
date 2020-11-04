/*
  subdomain.h: This file is part of Free Molecular Dynamics

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

#ifndef SUBDOMAIN_H
#define SUBDOMAIN_H

#include "config.h"
#include "types.h"
#include "array.h"

typedef struct _ParticleListItem ParticleListItem_t;

typedef ParticleListItem_t *cell_t;

typedef struct
{
    cell_t ***grid;                   /* must be equal to grid_array->data */
    fmd_array3D_t grid_array;
    int myrank;                       /* rank of the local process in MD_comm */
    int numprocs;                     /* number of processes in MD_comm */
    fmd_ituple_t is;                  /* position of subdomain in the subdomain grid */
    fmd_ituple_t rank_of_lower_subd;  /* rank of the neighbor processes */
    fmd_ituple_t rank_of_upper_subd;
    fmd_ituple_t ic_start;            /* width of margin, corresponds to the first */
                                      /* local index in the interior of the subdomain */
    fmd_ituple_t ic_stop;             /* first local index in the upper margin */
    fmd_utuple_t cell_num;            /* number of cells in subdomain, including margin */
    fmd_utuple_t cell_num_nonmarg;
    fmd_ituple_t ic_global_firstcell; /* global index of the first cell of the subdomain */
    unsigned NumberOfParticles;
    fmd_ituple_t r;                   /* r[d] = fmd_t.nc[d] % fmd_t.ns[d]; */
    fmd_ituple_t w;                   /* w[d] = fmd_t.nc[d] / fmd_t.ns[d]; */
} SubDomain_t;

typedef struct _fmd fmd_t;

void fmd_subd_init(fmd_t *md);
void fmd_subd_free(fmd_t *md);
void _fmd_convert_pos_to_subd_coord(fmd_t *md, fmd_rtuple_t pos, fmd_rtuple_t s);


#endif /* SUBDOMAIN_H */
