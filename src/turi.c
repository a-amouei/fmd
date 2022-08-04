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
#include "timer.h"
#include "events.h"
#include "general.h"
#include "h5.h"

typedef struct
{
    void *elms;            /* elements, an array */
    unsigned size;         /* number of elements */
} array_t;

static unsigned identify_tcell_processes_set(fmd_t *md, fmd_rtuple_t tcellh,
  fmd_ituple_t itc, int **pset);

/* the argument "psets" is a pointer to a single instance of array_t.
   the "elms" member of "psets" is an array of array_t's. */
static int find_this_pset_in_psets(int np, int pset[], array_t *psets)
{
    for (int i = 0; i < psets->size; i++)
    {
        array_t *arrt = &(((array_t *)(psets->elms))[i]);

        if (np == arrt->size)
        {
            int *pset2 = arrt->elms;
            int j;

            for (j=0; pset[j] == pset2[j] && j<np; j++) ;

            if (j==np) return i;
        }
    }

    return -1; /* it is not there! */
}

static int arrtcmp(const void *arrt1, const void *arrt2)
{
    const array_t *a1 = arrt1;
    const array_t *a2 = arrt2;
    int *pset1 = a1->elms;
    int *pset2 = a2->elms;

    if (pset1[0] < pset2[0])
        return -1;
    else
        if (pset1[0] > pset2[0])
            return 1;
        else
            if (a1->size > a2->size)
                return -1;
            else
                return 1;
}

/* the argument "psets" is a pointer to a single instance of array_t.
   the "elms" member of "psets" is an array of array_t's. */
static void make_psets_array(fmd_t *md, turi_t *t, array_t *psets)
{
    psets->size = 0;
    psets->elms = NULL;

    fmd_ituple_t itc;

    LOOP3D(itc, _fmd_ThreeZeros_int, t->tdims_global)
    {
        int *pset, np;

        np = identify_tcell_processes_set(md, t->tcellh, itc, &pset);

        if (np > 1)
        {
            int ip = find_this_pset_in_psets(np, pset, psets);

            if (ip == -1) /* if this is a new pset */
            {
                /* add this pset to psets */

                psets->elms = re_alloc(psets->elms, (psets->size+1) * sizeof(array_t));

                array_t *elms = psets->elms;

                elms[psets->size].elms = pset;
                elms[psets->size].size = np;

                psets->size++;
            }
            else
                free(pset);
        }
        else
            free(pset);
    }

    qsort(psets->elms, psets->size, sizeof(array_t), arrtcmp);
}

static fmd_bool_t is_table_row_available(array_t table[], int pset[], int sizepset, int irow)
{
    for (int i=0; i < sizepset; i++)
    {
        int icol = pset[i];

        if (table[icol].size >= irow+1)
            if ( ((int *)(table[icol].elms))[irow] != -1 )
                return FMD_FALSE;
    }

    return FMD_TRUE;
}

static void extend_table_column(array_t *col, unsigned newsize)
{
    col->elms = re_alloc(col->elms, newsize * sizeof(int));

    for (int i = col->size; i < newsize; i++)
        ((int *)(col->elms))[i] = -1;

    col->size = newsize;
}

static void insert_pset_in_table(array_t table[], int pset[], int sizepset, int psets_ind, int irow)
{
    for (int i=0; i < sizepset; i++)
    {
        int icol = pset[i];
        array_t *col = table + icol;

        if (col->size < irow+1) extend_table_column(col, irow+1);

        ((int *)(col->elms))[irow] = psets_ind;
    }
}

static void obtain_local_psets(fmd_t *md, turi_t *t, array_t *lpsets)
{
    /* create the (global) array of psets */

    array_t psets;

    make_psets_array(md, t, &psets);

    /* initialize the table */

    array_t *table;

    table = m_alloc(md->SubDomain.numprocs * sizeof(array_t));

    for (int i=0; i < md->SubDomain.numprocs; i++)
    {
        table[i].elms = NULL;
        table[i].size = 0;
    }

    /* fill the table */

    for (int psets_ind = 0; psets_ind < psets.size; psets_ind++)
    {
        array_t *pset = (array_t *)(psets.elms) + psets_ind;
        int *elpset = pset->elms;

        if (elpset[0] > md->SubDomain.myrank) break; /* because psets is sorted */

        for (int irow=0; ; irow++)
            if (is_table_row_available(table, elpset, pset->size, irow))
            {
                insert_pset_in_table(table, elpset, pset->size, psets_ind, irow);
                break;
            }
    }

    /* write output to lpsets */

    array_t *col = table + md->SubDomain.myrank;

    lpsets->size = 0;
    lpsets->elms = NULL;

    for (int i=0; i < col->size; i++)
    {
        int psets_ind = ((int *)(col->elms))[i];

        if (psets_ind != -1)
        {
            /* copy pset from psets to lpsets */

            lpsets->elms = re_alloc(lpsets->elms, (lpsets->size + 1) * sizeof(array_t));

            array_t *pset_src = (array_t *)(psets.elms) + psets_ind;
            array_t *pset_des = (array_t *)(lpsets->elms) + lpsets->size;

            pset_des->size = pset_src->size;
            pset_des->elms = m_alloc(pset_des->size * sizeof(int));

            memcpy(pset_des->elms, pset_src->elms, pset_des->size * sizeof(int));

            lpsets->size++;
        }
    }

    /* free memory allocated for "table" and "psets" */

    for (int i=0; i < psets.size; i++)
        free(((array_t *)(psets.elms))[i].elms);
    free(psets.elms);

    for (int i=0; i < md->SubDomain.numprocs; i++)
        free(table[i].elms);
    free(table);
}

static void call_field_update_event_handler(fmd_t *md, int field_index, int turi_index)
{
    fmd_event_params_field_update_t params;

    params.field = field_index;
    params.turi = turi_index;

    md->EventHandler(md, FMD_EVENT_FIELD_UPDATE, (fmd_event_params_t *)&params);
}

/* If the argument x is extremely close to an integer,
   returns the integer. Otherwise returns x unchanged. */
inline fmd_real_t impreal(fmd_real_t x)
{
    fmd_real_t e = fabs(x - round(x));

#ifdef FMD_FLOAT_CALCS
    return e < 1e-5 ? round(x) : x;
#else
    return e < 1e-10 ? round(x) : x;
#endif
}

static void gather_field_data_real(fmd_t *md, turi_t *t, field_t *f, fmd_array3D_t *out)
{
    out->data = NULL;
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* only "owners" participate in the collective communication below */
    if (owcomm->owned_tcells_num == 0) return;

    fmd_real_t ***lcd = (fmd_real_t ***)f->data.data;

    fmd_real_t *sendbuf = (fmd_real_t *)m_alloc(owcomm->owned_tcells_num * sizeof(*sendbuf));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        sendbuf[i] = lcd[itc[0]][itc[1]][itc[2]];
    }

    fmd_real_t *recvbuf;

    if (md->Is_MD_comm_root)
        recvbuf = (fmd_real_t *)m_alloc(t->tcells_global_num * sizeof(*recvbuf));

    MPI_Gatherv(sendbuf, owcomm->owned_tcells_num, FMD_MPI_REAL,
                recvbuf, owcomm->recvcounts, owcomm->displs,
                FMD_MPI_REAL, 0, owcomm->comm);

    free(sendbuf);

    if (md->Is_MD_comm_root)
    {
        _fmd_array_3d_create(t->tdims_global, sizeof(fmd_real_t), f->datatype, out);
        assert(out->data != NULL);
        /* TO-DO: handle memory error */

        for (int i=0; i < t->tcells_global_num; i++)
        {
            int *itc = owcomm->global_indexes[i];
            ((fmd_real_t ***)out->data)[itc[0]][itc[1]][itc[2]] = recvbuf[i];
        }

        free(recvbuf);
    }
}

static void gather_field_data_unsigned(fmd_t *md, turi_t *t, field_t *f, fmd_array3D_t *out)
{
    out->data = NULL;
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* only "owners" participate in the collective communication below */
    if (owcomm->owned_tcells_num == 0) return;

    unsigned ***lcd = (unsigned ***)f->data.data;

    unsigned *sendbuf = (unsigned *)m_alloc(owcomm->owned_tcells_num * sizeof(*sendbuf));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        sendbuf[i] = lcd[itc[0]][itc[1]][itc[2]];
    }

    unsigned *recvbuf;

    if (md->Is_MD_comm_root)
        recvbuf = (unsigned *)m_alloc(t->tcells_global_num * sizeof(*recvbuf));

    MPI_Gatherv(sendbuf, owcomm->owned_tcells_num, MPI_UNSIGNED,
                recvbuf, owcomm->recvcounts, owcomm->displs,
                MPI_UNSIGNED, 0, owcomm->comm);

    free(sendbuf);

    if (md->Is_MD_comm_root)
    {
        _fmd_array_3d_create(t->tdims_global, sizeof(unsigned), f->datatype, out);
        assert(out->data != NULL);
        /* TO-DO: handle memory error */

        for (int i=0; i < t->tcells_global_num; i++)
        {
            int *itc = owcomm->global_indexes[i];
            ((unsigned ***)out->data)[itc[0]][itc[1]][itc[2]] = recvbuf[i];
        }

        free(recvbuf);
    }
}

static void gather_field_data_rtuple(fmd_t *md, turi_t *t, field_t *f, fmd_array3D_t *out)
{
    out->data = NULL;
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* only "owners" participate in the collective communication below */
    if (owcomm->owned_tcells_num == 0) return;

    fmd_rtuple_t ***lcd = (fmd_rtuple_t ***)f->data.data;

    fmd_rtuple_t *sendbuf = (fmd_rtuple_t *)m_alloc(owcomm->owned_tcells_num * sizeof(*sendbuf));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        for (int d=0; d<3; d++)
            sendbuf[i][d] = lcd[itc[0]][itc[1]][itc[2]][d];
    }

    fmd_rtuple_t *recvbuf;

    if (md->Is_MD_comm_root)
        recvbuf = (fmd_rtuple_t *)m_alloc(t->tcells_global_num * sizeof(*recvbuf));

    MPI_Gatherv(sendbuf, owcomm->owned_tcells_num, md->mpi_types.mpi_rtuple,
                recvbuf, owcomm->recvcounts, owcomm->displs,
                md->mpi_types.mpi_rtuple, 0, owcomm->comm);

    free(sendbuf);

    if (md->Is_MD_comm_root)
    {
        _fmd_array_3d_create(t->tdims_global, sizeof(fmd_rtuple_t), f->datatype, out);
        assert(out->data != NULL);
        /* TO-DO: handle memory error */

        for (int i=0; i < t->tcells_global_num; i++)
        {
            int *itc = owcomm->global_indexes[i];

            for (int d=0; d<3; d++)
                ((fmd_rtuple_t ***)out->data)[itc[0]][itc[1]][itc[2]][d] = recvbuf[i][d];
        }

        free(recvbuf);
    }
}

/* This function receives the index of a turi-cell and returns the number
   of processes that share it. The ranks of those processes in MD_comm
   will be written in pset. */
static unsigned identify_tcell_processes_set(fmd_t *md, fmd_rtuple_t tcellh,
  fmd_ituple_t itc, int **pset)
{
    fmd_rtuple_t tc_edge_lo, tc_edge_hi;

    for (int d=0; d<3; d++)
    {
        tc_edge_lo[d] = itc[d] * tcellh[d];
        tc_edge_hi[d] = (itc[d] + 1) * tcellh[d];
    }

    fmd_rtuple_t slo, shi;
    _fmd_convert_pos_to_subd_coord(md, tc_edge_lo, slo);
    _fmd_convert_pos_to_subd_coord(md, tc_edge_hi, shi);

    fmd_ituple_t is_start, is_stop;
    unsigned np = 1; /* number of processes */

    for (int d=0; d<3; d++)
    {
        is_start[d] = (int)impreal(slo[d]);
        is_stop[d] = (int)ceil(impreal(shi[d]));
        np *= is_stop[d] - is_start[d];
    }

    *pset = (int *)m_alloc(np * sizeof(int));

    fmd_ituple_t is;
    int i=0;

    LOOP3D(is, is_start, is_stop)
        (*pset)[i++] = INDEX_FLAT(is, md->ns);

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

static void prepare_ownerscomm(fmd_t *md, turi_t *t)
{
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* create the communicator */

    MPI_Comm_split(md->MD_comm,
                   owcomm->owned_tcells_num == 0,
                   md->Is_MD_comm_root ? 0 : 1,
                   &owcomm->comm);

    /* prepare the remaining fields of t->ownerscomm */

    if (owcomm->owned_tcells_num == 0) return;

    MPI_Comm_size(owcomm->comm, &owcomm->commsize);

    if (md->Is_MD_comm_root)
        owcomm->recvcounts = (int *)m_alloc(owcomm->commsize * sizeof(int));

    MPI_Gather(&owcomm->owned_tcells_num, 1, MPI_INT, owcomm->recvcounts, 1, MPI_INT, 0, owcomm->comm);

    if (md->Is_MD_comm_root)
    {
        owcomm->global_indexes = (fmd_ituple_t *)m_alloc(t->tcells_global_num * sizeof(fmd_ituple_t));

        owcomm->displs = (int *)m_alloc(owcomm->commsize * sizeof(int));

        owcomm->displs[0] = 0;
        for (int i=1; i < owcomm->commsize; i++)
            owcomm->displs[i] = owcomm->displs[i-1] + owcomm->recvcounts[i-1];
    }

    fmd_ituple_t *lg = (fmd_ituple_t *)m_alloc(owcomm->owned_tcells_num * sizeof(fmd_ituple_t));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
        for (int d=0; d<3; d++)
            lg[i][d] = owcomm->owned_tcells[i][d] - t->itc_glob_to_loc[d];  /* convert local indexes to global indexes */

    MPI_Gatherv(lg, owcomm->owned_tcells_num, md->mpi_types.mpi_ituple, owcomm->global_indexes,
                owcomm->recvcounts, owcomm->displs, md->mpi_types.mpi_ituple, 0, owcomm->comm);

    free(lg);
}

static void prepare_turi_for_communication(fmd_t *md, turi_t *t)
{
    /* obtain local psets and create t->comms based on that */

    array_t lpsets;

    obtain_local_psets(md, t, &lpsets);

    if (lpsets.size > 0)
        t->comms = (turi_comm_t *)m_alloc(lpsets.size * sizeof(turi_comm_t));
    else
        t->comms = NULL;

    MPI_Group mdgroup, newgroup;

    MPI_Comm_group(md->MD_comm, &mdgroup);

    for (int i=0; i < lpsets.size; i++)
    {
        turi_comm_t *tcomm = t->comms + i;
        array_t *psetarrt = (array_t *)lpsets.elms + i;

        tcomm->commsize = psetarrt->size;
        tcomm->pset = psetarrt->elms;
        tcomm->num_tcells = 0;
        tcomm->itcs = NULL;

        /* create MPI communicator */
        MPI_Group_incl(mdgroup, tcomm->commsize, tcomm->pset, &newgroup);
        int res = MPI_Comm_create_group(md->MD_comm, newgroup, 0, &tcomm->comm);
        assert(res == MPI_SUCCESS);
        MPI_Group_free(&newgroup);
    }

    t->comms_num = lpsets.size;

    free(lpsets.elms);

    MPI_Group_free(&mdgroup);

    /* find "owned" turi-cells and associate each individual turi-cell with one of the turi_comm_t's if the
       turi-cell is shared with other subdomains */

    t->ownerscomm.owned_tcells = NULL;
    t->ownerscomm.owned_tcells_num = 0;

    fmd_ituple_t itc;

    LOOP3D(itc, t->itc_start_global, t->itc_stop_global)  /* a loop on all turi-cells in current subdomain */
    {
        int *pset, np;

        np = identify_tcell_processes_set(md, t->tcellh, itc, &pset);

        /* is this turi-cell owned by current subdomain? (do owned_tcells and owned_tcells_num need an update?) */

        if (pset[0] == md->SubDomain.myrank)
        {
            t->ownerscomm.owned_tcells = (fmd_ituple_t *)re_alloc(t->ownerscomm.owned_tcells,
                                          (t->ownerscomm.owned_tcells_num+1) * sizeof(*t->ownerscomm.owned_tcells));

            for (int d=0; d<3; d++)
                t->ownerscomm.owned_tcells[t->ownerscomm.owned_tcells_num][d] = itc[d] + t->itc_glob_to_loc[d];

            t->ownerscomm.owned_tcells_num++;
        }

        /* if this turi-cell is shared with other subdomains, associate it with one of the turi_comm_t's */

        if (np > 1)
        {
            int icomm = find_this_pset_in_comms(np, pset, t->comms_num, t->comms);

            assert(icomm > -1);

            turi_comm_t *tcomm = t->comms + icomm;

            /* now, add the local index of the current turi-cell to "itcs" array */

            tcomm->itcs = (fmd_ituple_t *)re_alloc(tcomm->itcs, (tcomm->num_tcells+1) * sizeof(fmd_ituple_t));

            for (int d=0; d<3; d++)
                tcomm->itcs[tcomm->num_tcells][d] = itc[d] + t->itc_glob_to_loc[d];

            tcomm->num_tcells++;

        }

        free(pset);
    }

    /* complete preparation of ownerscomm */

    prepare_ownerscomm(md, t);

    /* find num_tcells_max */

    t->num_tcells_max = 0;
    for (int i=0; i < t->comms_num; i++)
        if (t->comms[i].num_tcells > t->num_tcells_max) t->num_tcells_max = t->comms[i].num_tcells;
}

/* find in each direction the ranks of the upper and lower subdomains which are "owners" of
   turi-cells, and initialize rank_of_lower_owner and rank_of_upper_owner accordingly. */
static void init_rank_of_lower_upper_owner(fmd_t *md, turi_t *t)
{
    fmd_ituple_t is_tempo;

    for (int d=0; d<DIM; d++)                 /* set is_tempo[] to current subdomain */
        is_tempo[d] = md->SubDomain.is[d];

    for (int d=0; d<DIM; d++)
    {
        fmd_real_t pos;

        /* initialize rank_of_upper_owner[d] */

        pos = (t->itc_stop_global[d] % t->tdims_global[d]) * t->tcellh[d];
        is_tempo[d] = (int)impreal( _fmd_convert_pos_to_subd_coord_1D(md, pos, d) );
        t->rank_of_upper_owner[d] = INDEX_FLAT(is_tempo, md->ns);

        if ( (t->rank_of_upper_owner[d] < md->SubDomain.myrank && !md->PBC[d]) ||
             (t->rank_of_upper_owner[d] == md->SubDomain.myrank) )
            t->rank_of_upper_owner[d] = MPI_PROC_NULL;

        /* initialize rank_of_lower_owner[d] */

        int ifst = t->itc_start_owned[d] - t->itc_glob_to_loc[d]; /* convert to global index */
        pos = ((ifst - 1 + t->tdims_global[d]) % t->tdims_global[d]) * t->tcellh[d];
        is_tempo[d] = (int)impreal( _fmd_convert_pos_to_subd_coord_1D(md, pos, d) );
        t->rank_of_lower_owner[d] = INDEX_FLAT(is_tempo, md->ns);

        if ( (t->rank_of_lower_owner[d] > md->SubDomain.myrank && !md->PBC[d]) ||
             (t->rank_of_lower_owner[d] == md->SubDomain.myrank) )
            t->rank_of_lower_owner[d] = MPI_PROC_NULL;

        /* initialize has_upper_lower_procs[d] */

        t->has_upper_lower_owner_procs[d] = (t->rank_of_upper_owner[d] != MPI_PROC_NULL) ||
                                            (t->rank_of_lower_owner[d] != MPI_PROC_NULL);

        /* set is_tempo[] back again to current subdomain so as to prepare for the next iteration */

        is_tempo[d] = md->SubDomain.is[d];
    }
}

fmd_handle_t fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz, fmd_real_t starttime, fmd_real_t stoptime)
{
    if (md->SubDomain.grid == NULL) fmd_subd_init(md);

    int ti = md->turies_num;

    md->turies = (turi_t *)re_alloc(md->turies, (ti+1) * sizeof(turi_t));

    turi_t *t = &md->turies[ti];

    t->turi_index = ti;

    t->tdims_global[0] = dimx;
    t->tdims_global[1] = dimy;
    t->tdims_global[2] = dimz;

    t->starttime = starttime;
    t->stoptime = stoptime;

    t->tcells_global_num = dimx * dimy * dimz;

    for (int d=0; d<DIM; d++)
    {
        t->itc_start[d] = (cat == FMD_TURI_TTM) ? 1 : 0;

        t->tcellh[d] = md->l[d] / t->tdims_global[d];

        fmd_real_t xlo = md->SubDomain.ic_global_firstcell[d] * md->cellh[d];
        fmd_real_t xlodiv = impreal(xlo / t->tcellh[d]);
        t->itc_start_global[d] = (int)xlodiv;

        /* if the fractional part of xlodiv is zero, then t->itc_start_owned[d] = t->itc_start[d] */
        t->itc_start_owned[d] = t->itc_start[d] + ((fmd_real_t)t->itc_start_global[d] == xlodiv ? 0 : 1);

        fmd_real_t xhi = xlo + md->SubDomain.cell_num_nonmarg[d] * md->cellh[d];
        t->itc_stop_global[d] = (int)ceil(impreal(xhi / t->tcellh[d]));

        t->tdims_local_nonmarg[d] = t->itc_stop_global[d] - t->itc_start_global[d];

        t->itc_stop[d] = t->itc_start[d] + t->tdims_local_nonmarg[d];

        t->tdims_local[d] = t->itc_start[d] + t->itc_stop[d];

        t->itc_glob_to_loc[d] = -t->itc_start_global[d] + t->itc_start[d];
    }

    t->tcell_volume = t->tcellh[0] * t->tcellh[1] * t->tcellh[2];

    t->cat = cat;

    t->fields = NULL;
    t->fields_num = 0;

    switch (cat)
    {
        case FMD_TURI_CUSTOM:
            break;

        case FMD_TURI_TTM:
            break;

        default:
            assert(0);  /* TO-DO */
    }

    init_rank_of_lower_upper_owner(md, t);

    md->turies_num++;

    prepare_turi_for_communication(md, t);

    return ti;
}

static void set_field_data_el_size_and_type(field_t *f)
{
    if (f->cat == FMD_FIELD_MASS || f->cat == FMD_FIELD_TEMPERATURE || f->cat == FMD_FIELD_NUMBER_DENSITY)
    {
        f->data_el_size = sizeof(fmd_real_t);
        f->datatype = DATATYPE_REAL;
    }
    else if (f->cat == FMD_FIELD_VCM)
    {
        f->data_el_size = sizeof(fmd_rtuple_t);
        f->datatype = DATATYPE_RTUPLE;
    }
    else if (f->cat == FMD_FIELD_NUMBER)
    {
        f->data_el_size = sizeof(unsigned);
        f->datatype = DATATYPE_UNSIGNED;
    }
    else
        assert(0);  /* TO-DO */
}

inline static void prepare_buf1_real(fmd_real_t *buf1, turi_comm_t *tcm, field_t *f)
{
    for (int j=0; j < tcm->num_tcells; j++)
    {
        int *itc = tcm->itcs[j];

        buf1[j] = ((fmd_real_t ***)f->data.data)[itc[0]][itc[1]][itc[2]];
    }
}

inline static void prepare_buf1_rtuple(fmd_rtuple_t *buf1, turi_comm_t *tcm, field_t *f)
{
    for (int j=0; j < tcm->num_tcells; j++)
    {
        int *itc = tcm->itcs[j];

        for (int d=0; d<3; d++)
            buf1[j][d] = ((fmd_rtuple_t ***)f->data.data)[itc[0]][itc[1]][itc[2]][d];
    }
}

inline static void prepare_buf1_unsigned(unsigned *buf1, turi_comm_t *tcm, field_t *f)
{
    for (int j=0; j < tcm->num_tcells; j++)
    {
        int *itc = tcm->itcs[j];

        buf1[j] = ((unsigned ***)f->data.data)[itc[0]][itc[1]][itc[2]];
    }
}

static void perform_field_comm_rtuple(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    /* allocate memory for buf1 and buf2 */
    fmd_rtuple_t *buf1 = m_alloc(t->num_tcells_max * sizeof(*buf1));
    fmd_rtuple_t *buf2 = m_alloc(t->num_tcells_max * sizeof(*buf2));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buf1 (send-buffer) */

        prepare_buf1_rtuple(buf1, tcm, f);

        /* do communication */

        if (allreduce)
            MPI_Allreduce(buf1, buf2, 3*tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, tcm->comm);
        else
            MPI_Reduce(buf1, buf2, 3*tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, 0, tcm->comm);

        /* copy from buf2 to data array */

        if (allreduce || tcm->pset[0] == md->SubDomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                for (int d=0; d<3; d++)
                    ((fmd_rtuple_t ***)f->data.data)[itc[0]][itc[1]][itc[2]][d] = ((fmd_rtuple_t *)buf2)[j][d];
            }
        }
    }

    free(buf1);
    free(buf2);
}

static void perform_field_comm_real(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    /* allocate memory for buf1 and buf2 */
    fmd_real_t *buf1 = m_alloc(t->num_tcells_max * sizeof(*buf1));
    fmd_real_t *buf2 = m_alloc(t->num_tcells_max * sizeof(*buf2));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buf1 (send-buffer) */

        prepare_buf1_real(buf1, tcm, f);

        /* do communication */

        if (allreduce)
            MPI_Allreduce(buf1, buf2, tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, tcm->comm);
        else
            MPI_Reduce(buf1, buf2, tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, 0, tcm->comm);

        /* copy from buf2 to data array */

        if (allreduce || tcm->pset[0] == md->SubDomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                ((fmd_real_t ***)f->data.data)[itc[0]][itc[1]][itc[2]] = ((fmd_real_t *)buf2)[j];
            }
        }
    }

    free(buf1);
    free(buf2);
}

static void perform_field_comm_unsigned(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    /* allocate memory for buf1 and buf2 */
    unsigned *buf1 = m_alloc(t->num_tcells_max * sizeof(*buf1));
    unsigned *buf2 = m_alloc(t->num_tcells_max * sizeof(*buf2));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buf1 (send-buffer) */

        prepare_buf1_unsigned(buf1, tcm, f);

        /* do communication */

        if (allreduce)
            MPI_Allreduce(buf1, buf2, tcm->num_tcells, MPI_UNSIGNED, MPI_SUM, tcm->comm);
        else
            MPI_Reduce(buf1, buf2, tcm->num_tcells, MPI_UNSIGNED, MPI_SUM, 0, tcm->comm);

        /* copy from buf2 to data array */

        if (allreduce || tcm->pset[0] == md->SubDomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                ((unsigned ***)f->data.data)[itc[0]][itc[1]][itc[2]] = ((unsigned *)buf2)[j];
            }
        }

    }

    free(buf1);
    free(buf2);
}

static void update_field_number(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    unsigned ***num = (unsigned ***)f->data.data;

    /* clean data of number field (initialize with zeros) */
    _fmd_array_3d_unsigned_clean(num, t->tdims_local);

    fmd_ituple_t ic, itc;
    cell_t *cell;
    int i;
    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {
            /* find local index of turi-cell from particle's position */
            for (int d=0; d<3; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* count atoms */
            num[itc[0]][itc[1]][itc[2]]++;
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_comm_unsigned(md, f, t, allreduce);

    f->timestep = md->time_iteration; /* mark as updated */

    if (md->EventHandler != NULL) call_field_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_number_density(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    field_t *fdep = &t->fields[f->dependcs[0]];
    /* update the dependency field if not already updated */
    if (fdep->timestep != md->time_iteration) update_field_number(md, fdep, t, allreduce || fdep->allreduce_now);

    unsigned ***num = (unsigned ***)fdep->data.data;
    fmd_real_t ***nd = (fmd_real_t ***)f->data.data;

    if (allreduce)
    {
        fmd_ituple_t itc;

        LOOP3D(itc, t->itc_start, t->itc_stop)
            nd[itc[0]][itc[1]][itc[2]] = num[itc[0]][itc[1]][itc[2]] / t->tcell_volume;
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            nd[itc[0]][itc[1]][itc[2]] = num[itc[0]][itc[1]][itc[2]] / t->tcell_volume;
        }
    }

    f->timestep = md->time_iteration; /* mark as updated */

    if (md->EventHandler != NULL) call_field_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_mass(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    fmd_real_t ***mass = (fmd_real_t ***)f->data.data;

    /* clean data of mass field (initialize with zeros) */
    _fmd_array_3d_real_clean(mass, t->tdims_local);

    fmd_ituple_t ic, itc;
    cell_t *cell;
    int i;
    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {
            /* find index of turi-cell from particle's position */
            for (int d=0; d<3; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* update mass for turi-cell with index of itc */
            mass[itc[0]][itc[1]][itc[2]] += md->potsys.atomkinds[cell->atomkind[i]].mass;
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_comm_real(md, f, t, allreduce);

    f->timestep = md->time_iteration; /* mark as updated */

    if (md->EventHandler != NULL) call_field_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_vcm_only(fmd_t *md, field_t *fvcm, field_t *fmass, turi_t *t, fmd_bool_t allreduce)
{
    fmd_ituple_t ic, itc;

    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)fvcm->data.data;

    /* clean data of vcm field (initialize with zeros) */
    _fmd_array_3d_rtuple_clean(vcm, t->tdims_local);

    cell_t *cell;
    int i;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {

            /* find index of turi-cell from particle's position */
            for (int d=0; d<3; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* calculate vcm and mass (first phase) */

            fmd_real_t m = md->potsys.atomkinds[cell->atomkind[i]].mass;

            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] += m * VEL(cell, i, d);
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_comm_rtuple(md, fvcm, t, allreduce);

    /* calculate vcm (last phase) */

    fmd_real_t ***mass = (fmd_real_t ***)fmass->data.data;

    if (allreduce)
    {
        LOOP3D(itc, t->itc_start, t->itc_stop)
            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] /= mass[itc[0]][itc[1]][itc[2]];
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] /= mass[itc[0]][itc[1]][itc[2]];
        }
    }

    fvcm->timestep = md->time_iteration; /* mark as updated */
    if (md->EventHandler != NULL) call_field_update_event_handler(md, fvcm->field_index, t->turi_index);
}

static void update_field_vcm_and_mass(fmd_t *md, field_t *fvcm, field_t *fmass, turi_t *t, fmd_bool_t allreduce)
{
    fmd_ituple_t ic, itc;

    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)fvcm->data.data;
    fmd_real_t ***mass = (fmd_real_t ***)fmass->data.data;

    /* clean data of mass and vcm fields (initialize with zeros) */
    LOOP3D(itc, t->itc_start, t->itc_stop)
    {
        mass[itc[0]][itc[1]][itc[2]] = 0.0;
        for (int d=0; d<3; d++)
            vcm[itc[0]][itc[1]][itc[2]][d] = 0.0;
    }

    cell_t *cell;
    int i;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {

            /* find index of turi-cell from particle's position */
            for (int d=0; d<3; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* calculate vcm and mass (first phase) */

            fmd_real_t m = md->potsys.atomkinds[cell->atomkind[i]].mass;

            mass[itc[0]][itc[1]][itc[2]] += m;
            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] += m * VEL(cell, i, d);
        }

    /* do communications */
    if (t->comms_num > 0)
    {
        perform_field_comm_real(md, fmass, t, allreduce || fmass->allreduce_now);
        perform_field_comm_rtuple(md, fvcm, t, allreduce);
    }

    /* calculate vcm (last phase) */
    if (allreduce)
    {
        LOOP3D(itc, t->itc_start, t->itc_stop)
            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] /= mass[itc[0]][itc[1]][itc[2]];
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] /= mass[itc[0]][itc[1]][itc[2]];
        }
    }

    fmass->timestep = md->time_iteration; /* mark as updated */
    if (md->EventHandler != NULL) call_field_update_event_handler(md, fmass->field_index, t->turi_index);
    fvcm->timestep = md->time_iteration; /* mark as updated */
    if (md->EventHandler != NULL) call_field_update_event_handler(md, fvcm->field_index, t->turi_index);
}

static void update_field_vcm(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    field_t *fmass = &t->fields[f->dependcs[0]];

    if (fmass->timestep != md->time_iteration)  /* if fmass is not already updated */
        update_field_vcm_and_mass(md, f, fmass, t, allreduce);
    else
        update_field_vcm_only(md, f, fmass, t, allreduce);
}

#define ADD_interval_TO_ARRAY(array, num, interval)                            \
    do                                                                         \
    {                                                                          \
        (array) = (fmd_real_t *)re_alloc(array, ((num)+1)*sizeof(fmd_real_t)); \
        (array)[(num)++] = (interval);                                         \
    } while(0)

static void field_intervals_add(field_t *f, fmd_real_t interval, fmd_bool_t allreduce)
{
    unsigned j;

    /* see if the interval already exists in "intervals_allreduce" array */
    for (j=0; j < f->intervals_allreduce_num; j++)
        if (f->intervals_allreduce[j] == interval) break;

    if (j == f->intervals_allreduce_num) /* no, doesn't exist there */
    {
        /* see if the interval already exists in "intervals" array */
        for (j=0; j < f->intervals_num; j++)
            if (f->intervals[j] == interval) break;

        if (j == f->intervals_num) /* doesn't exist in "intervals" either */
        {
            if (allreduce) /* add it to "intervals_allreduce" */
                ADD_interval_TO_ARRAY(f->intervals_allreduce, f->intervals_allreduce_num, interval);
            else /* add it to "intervals" */
                ADD_interval_TO_ARRAY(f->intervals, f->intervals_num, interval);
        }
        else /* it exists in "intervals" array */
        {
            if (allreduce)
            {
                /* remove it from "intervals" */
                f->intervals[j] = --f->intervals_num;

                /* add it to "intervals_allreduce" */
                ADD_interval_TO_ARRAY(f->intervals_allreduce, f->intervals_allreduce_num, interval);
            }
        }
    }
}

static fmd_handle_t field_add(fmd_t *md, fmd_handle_t turi, fmd_field_t cat, fmd_real_t interval, fmd_bool_t allreduce)
{
    int i;
    field_t *f;
    unsigned dep1;
    turi_t *t = &md->turies[turi];

    /* add dependency fields */
    switch (cat)
    {
        case FMD_FIELD_TEMPERATURE:
            dep1 = field_add(md, turi, FMD_FIELD_VCM, interval, FMD_TRUE);
            break;

        case FMD_FIELD_VCM:
            dep1 = field_add(md, turi, FMD_FIELD_MASS, interval, allreduce);
            break;

        case FMD_FIELD_NUMBER_DENSITY:
            dep1 = field_add(md, turi, FMD_FIELD_NUMBER, interval, allreduce);
            break;

        default:
            ;
    }

    /* check if this field is already added */
    for (i=0; i < t->fields_num; i++)
        if (t->fields[i].cat == cat) break;

    /* if not, add it */
    if (i == t->fields_num)
    {
        /* add the field */

        t->fields = (field_t *)re_alloc(t->fields, (t->fields_num+1) * sizeof(field_t));

        f = &t->fields[t->fields_num];
        f->field_index = t->fields_num;
        f->cat = cat;
        f->timestep = -1;
        f->intervals_allreduce = NULL;
        f->intervals_allreduce_num = 0;
        f->intervals = NULL;
        f->intervals_num = 0;
        t->fields_num++;

        if (cat == FMD_FIELD_TEMPERATURE || cat == FMD_FIELD_VCM || cat == FMD_FIELD_NUMBER_DENSITY)
        {
            f->dependcs_num = 1;
            f->dependcs = (unsigned *)m_alloc(sizeof(unsigned));
            f->dependcs[0] = dep1;
        }
        else
            f->dependcs_num = 0;

        set_field_data_el_size_and_type(f);

        /* allocate space for field data */
        _fmd_array_3d_create(t->tdims_local, f->data_el_size, f->datatype, &f->data);
        assert(f->data.data != NULL);
    }
    else
        f = &t->fields[i];

    /* add the interval if doesn't already exist in "intervals" or "intervals_allreduce" arrays */
    field_intervals_add(f, interval, allreduce);

    return i;
}

fmd_handle_t fmd_field_add(fmd_t *md, fmd_handle_t turi, fmd_field_t cat, fmd_real_t interval)
{
    return field_add(md, turi, cat, interval, FMD_FALSE);
}

static void update_field(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    switch (f->cat)
    {
        case FMD_FIELD_NUMBER:
            update_field_number(md, f, t, allreduce);
            break;

        case FMD_FIELD_NUMBER_DENSITY:
            update_field_number_density(md, f, t, allreduce);
            break;

        case FMD_FIELD_MASS:
            update_field_mass(md, f, t, allreduce);
            break;

        case FMD_FIELD_VCM:
            update_field_vcm(md, f, t, allreduce);
            break;
    }
}

/* this subroutine updates the variable allreduce_now in all fields */
static inline void allreduce_now_update(fmd_t *md, turi_t *t)
{
    for (int fi = 0; fi < t->fields_num; fi++)
    {
        field_t *f = &t->fields[fi];

        f->allreduce_now = FMD_FALSE;

        for (int j = 0; j < f->intervals_allreduce_num; j++)
            if (_fmd_timer_is_its_time(md->MD_time, md->delta_t/2.0, t->starttime, f->intervals_allreduce[j]))
            {
                f->allreduce_now = FMD_TRUE;
                break;
            }
    }
}

void _fmd_turies_update(fmd_t *md)
{
    for (int ti=0; ti < md->turies_num; ti++)
    {
        turi_t *t = &md->turies[ti];

        if (md->MD_time >= t->starttime && !(md->MD_time > t->stoptime && t->stoptime >= t->starttime))
        {
            allreduce_now_update(md, t);

            /* start from fields with higher indexes */
            for (int fi = t->fields_num-1; fi >= 0; fi--)
            {
                field_t *f = &t->fields[fi];

                if (f->timestep == md->time_iteration) continue; /* see if the field is already updated */

                if (f->allreduce_now)
                {
                    update_field(md, f, t, FMD_TRUE);
                    break;
                }
                else
                {
                    for (int j=0; j < f->intervals_num; j++)
                        if (_fmd_timer_is_its_time(md->MD_time, md->delta_t/2.0, t->starttime, f->intervals[j]))
                        {
                            update_field(md, f, t, FMD_FALSE);
                            break;
                        }
                }
            }
        }
    }
}

static void convert_inmass_to_outmass(fmd_array3D_t *ar)
{
    fmd_utriple_t iv;

    LOOP3D(iv, _fmd_ThreeZeros_int, ar->dims)
        ((fmd_real_t ***)(ar->data))[iv[0]][iv[1]][iv[2]] *= MD_MASS_UNIT;
}

void fmd_field_save_as_hdf5(fmd_t *md, fmd_handle_t turi, fmd_handle_t field, fmd_string_t path)
{
    turi_t *t = &md->turies[turi];
    field_t *f = &t->fields[field];

    fmd_array3D_t data;

    switch(f->cat)
    {
        case FMD_FIELD_VCM:
            gather_field_data_rtuple(md, t, f, &data);
            if (md->Is_MD_comm_root)
            {
                _fmd_h5_save_tuple_field_float(md, "vcm", t, path, &data);
                _fmd_array_3d_free(&data);
            }
            break;

        case FMD_FIELD_NUMBER_DENSITY:
            gather_field_data_real(md, t, f, &data);
            if (md->Is_MD_comm_root)
            {
                _fmd_h5_save_scalar_field_float(md, "number-density", t, path, &data);
                _fmd_array_3d_free(&data);
            }
            break;

        case FMD_FIELD_MASS:
            gather_field_data_real(md, t, f, &data);
            if (md->Is_MD_comm_root)
            {
                convert_inmass_to_outmass(&data);
                _fmd_h5_save_scalar_field_float(md, "mass", t, path, &data);
                _fmd_array_3d_free(&data);
            }
            break;

        default:
            assert(0);
    }
}

static void field_free(field_t *f)
{
    if (f->dependcs_num > 0) free(f->dependcs);
    free(f->intervals_allreduce);
    free(f->intervals);
    _fmd_array_3d_free(&f->data);
}

static void turi_comm_free(turi_comm_t *tcm)
{
    if (tcm->commsize > 1)
    {
        int res = MPI_Comm_free(&tcm->comm);
        assert(res == MPI_SUCCESS);
    }
    free(tcm->pset);
    free(tcm->itcs);
}

static void turi_ownerscomm_free(fmd_t *md, turi_ownerscomm_t *towcm)
{
    free(towcm->owned_tcells);

    int res = MPI_Comm_free(&towcm->comm);
    assert(res == MPI_SUCCESS);

    if (md->Is_MD_comm_root)
    {
        free(towcm->global_indexes);
        free(towcm->recvcounts);
        free(towcm->displs);
    }
}

static void turi_free(fmd_t *md, turi_t *t)
{
    for (unsigned u=0; u < t->fields_num; u++)
        field_free(&t->fields[u]);
    free(t->fields);

    for (unsigned u=0; u < t->comms_num; u++)
        turi_comm_free(&t->comms[u]);
    free(t->comms);

    turi_ownerscomm_free(md, &t->ownerscomm);
}

void fmd_turi_free(fmd_t *md)
{
    for (unsigned u=0; u < md->turies_num; u++)
        turi_free(md, &md->turies[u]);
    free(md->turies);

    md->turies = NULL;
    md->turies_num = 0;
}
