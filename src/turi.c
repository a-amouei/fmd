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

static int find_this_pset_in_comms(int np, int *pset, int comms_num, turi_comm_t *comms)
{
    for (int i=0; i<comms_num; i++)
        if (np == comms[i].commsize)
        {
            int *pset2 = comms[i].pset;
            int j;

            for (j=0; pset[j] == pset2[j] && j<np; j++) ;

            if (j==np) return i;
        }

    return -1; /* it is not there! */
}

static void prepare_for_communication(fmd_t *md, turi_t *t)
{
    t->comms_num = 0;
    t->comms = NULL;

    MPI_Group worldgroup, newgroup;

    MPI_Comm_group(MPI_COMM_WORLD, &worldgroup);

    int itc[3];

    ITERATE(itc, t->tcell_start, t->tcell_stop)
    {
        turi_comm_t *tcomm;
        int *pset, np;

        np = identify_tcell_processes_set(md, t->tcellh, itc, &pset);
        int icomm = find_this_pset_in_comms(np, pset, t->comms_num, t->comms);

        if (icomm == -1) /* if this is a new pset */
        {
            /* add it to t->comms */

            t->comms = (turi_comm_t *)realloc(t->comms, (t->comms_num+1) * sizeof(turi_comm_t));
            /* TO-DO: handle memory error */
            assert(t->comms != NULL);

            tcomm = &t->comms[t->comms_num++];
            tcomm->commsize = np;
            tcomm->pset = pset;
            tcomm->num_tcells = 0;
            tcomm->itcs = NULL;

            if (np > 1) /* create MPI communicator */
            {
                MPI_Group_incl(worldgroup, np, pset, &newgroup);
                int res = MPI_Comm_create_group(md->MD_comm, newgroup, 0, &tcomm->comm);
                assert(res == MPI_SUCCESS);
                MPI_Group_free(&newgroup);
            }
        }
        else
        {
            free(pset);
            tcomm = &t->comms[icomm];
        }

        /* now, add the local index of the current turi-cell to "itcs" array */

        tcomm->itcs = (index_t *)realloc(tcomm->itcs, (tcomm->num_tcells+1) * sizeof(index_t));
        /* TO-DO: handle memory error */
        assert(tcomm->itcs != NULL);

        for (int d=0; d<3; d++)
            tcomm->itcs[tcomm->num_tcells][d] = itc[d] - t->tcell_start[d];

        tcomm->num_tcells++;
    }

    MPI_Group_free(&worldgroup);
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

    prepare_for_communication(md, t);

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
