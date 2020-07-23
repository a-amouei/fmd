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

typedef struct _ParticleListItem ParticleListItem_t;

typedef ParticleListItem_t *cell_t;

typedef struct
{
    cell_t ***grid;             // where the particles lie
    int myrank;                 // rank of the local process in MD_comm
    int numprocs;               // number of processes in MD_comm
    int is[3];                  // position of subdomain in the subdomain grid
    int rank_of_lower_subd[3];  // rank of the neighbor processes
    int rank_of_upper_subd[3];
    int ic_start[3];            // width of margin, corresponds to the first
                                // local index in the interior of the subdomain
    int ic_stop[3];             // first local index in the upper margin
    int cell_num[3];            // number of cells in subdomain, including margin
    unsigned cell_num_nonmarg[3];
    int ic_global_firstcell[3]; // global index of the first cell of the subdomain
    unsigned NumberOfParticles;
    int r[3];                   /* r[d] = fmd_t.nc[d] % fmd_t.ns[d]; */
    int w[3];                   /* w[d] = fmd_t.nc[d] / fmd_t.ns[d]; */
} SubDomain_t;

typedef struct _fmd fmd_t;

void fmd_subd_init(fmd_t *md);
void fmd_subd_free(fmd_t *md);
void _fmd_convert_pos_to_subd_coord(fmd_t *md, const double pos[3], double s[3]);


#endif /* SUBDOMAIN_H */
