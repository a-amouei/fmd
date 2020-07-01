/*
  turi.c: This file is part of Free Molecular Dynamics

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

#include "turi.h"
#include "base.h"
#include "subdomain.h"

/* This function receives the index of a turi-cell and returns the number
   of processes that share it. The ranks of those processes in MD_comm
   will be written in pset. */
static unsigned identify_tcell_processes_set(fmd_t *md, double tcellh[3],
  const int itc[3], int *pset[])
{
    double tc_edge_lo[3], tc_edge_hi[3];

    for (int d=0; d<3; d++)
    {
        tc_edge_lo[d] = itc[d] * tcellh[d];
        tc_edge_hi[d] = (itc[d] + 1) * tcellh[d];
    }

    double slo[3], shi[3];
    _fmd_convert_pos_to_subd_coord(md, tc_edge_lo, slo);
    _fmd_convert_pos_to_subd_coord(md, tc_edge_hi, shi);

    int is_start[3], is_stop[3];
    unsigned np = 1; /* number of processes */

    for (int d=0; d<3; d++)
    {
        is_start[d] = (int)slo[d];
        is_stop[d] = (int)ceil(shi[d]);
        np *= is_stop[d] - is_start[d];
    }

    *pset = (int *)malloc(np * sizeof(int));
    /* TO-DO: handle memory error */
    assert(*pset != NULL);

    int is[3];
    int i=0;

    ITERATE(is, is_start, is_stop)
        (*pset)[i++] = INDEX(is, md->ns);

    return np;
}

unsigned fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz)
{
    if (md->SubDomain.grid == NULL) fmd_subd_init(md);

    int ti = md->turies_num;

    md->turies = (turi_t *)realloc(md->turies, (ti+1) * sizeof(turi_t));
    /* TO-DO: handle memory error */
    assert(md->turies != NULL);

    turi_t *t = &md->turies[ti];

    t->tdims_global[0] = dimx;
    t->tdims_global[1] = dimy;
    t->tdims_global[2] = dimz;

    for (int d=0; d<3; d++)
    {
        t->tcellh[d] = md->l[d] / t->tdims_global[d];

        double xlo = md->SubDomain.ic_global_firstcell[d] * md->cellh[d];
        t->tcell_start[d] = (int)(xlo / t->tcellh[d]);
        double xhi = xlo + md->SubDomain.cell_num_nonmarg[d] * md->cellh[d];
        t->tcell_stop[d] = (int)ceil(xhi / t->tcellh[d]);

        t->tdims[d] = t->tcell_stop[d] - t->tcell_start[d];
    }

    t->cat = cat;
    switch (cat)
    {
        case FMD_TURI_CUSTOM:
            t->fields = NULL;
            t->fields_num = 0;
            break;
    }

    md->turies_num++;

    /* only for test purpose */
    /*
    for (int i=0; i < md->ns[0]; i++)
        for (int j=0; j < md->ns[1]; j++)
            for (int k=0; k < md->ns[2]; k++)
            {
                MPI_Barrier(md->MD_comm);
                if (md->SubDomain.is[0] == i &&
                    md->SubDomain.is[1] == j &&
                    md->SubDomain.is[2] == k)
                {
                    printf("subdomain(%d, %d, %d)\n", i, j, k);
                    printf("tcell_start = {%d, %d, %d}\n", t->tcell_start[0], t->tcell_start[1], t->tcell_start[2]);
                    printf("tcell_stop = {%d, %d, %d}\n", t->tcell_stop[0], t->tcell_stop[1], t->tcell_stop[2]);
                    printf("tdims = {%d, %d, %d}\n\n", t->tdims[0], t->tdims[1], t->tdims[2]);
                }
            }
    */

    return ti;
}
