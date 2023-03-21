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

#include <string.h>
#include <tgmath.h>
#include "turi.h"
#include "fmd-private.h"
#include "misc.h"
#include "subdomain.h"
#include "timer.h"
#include "events.h"
#include "general.h"
#include "h5.h"
#include "ttm.h"

typedef struct
{
    void *elms;            /* elements, an array */
    unsigned size;         /* number of elements */
} array_t;

fmd_bool_t _is_time_within_turi_start_stop_times(fmd_t *md, turi_t *t)
{
    return (md->time >= t->starttime - md->timestep/2.0) &&
           ((md->time < t->stoptime + md->timestep/2.0) || t->stoptime < t->starttime);
}

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

    table = m_alloc(md->Subdomain.numprocs * sizeof(array_t));

    for (int i=0; i < md->Subdomain.numprocs; i++)
    {
        table[i].elms = NULL;
        table[i].size = 0;
    }

    /* fill the table */

    for (int psets_ind = 0; psets_ind < psets.size; psets_ind++)
    {
        array_t *pset = (array_t *)(psets.elms) + psets_ind;
        int *elpset = pset->elms;

        if (elpset[0] > md->Subdomain.myrank) break; /* because psets is sorted */

        for (int irow=0; ; irow++)
            if (is_table_row_available(table, elpset, pset->size, irow))
            {
                insert_pset_in_table(table, elpset, pset->size, psets_ind, irow);
                break;
            }
    }

    /* write output to lpsets */

    array_t *col = table + md->Subdomain.myrank;

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

    for (int i=0; i < md->Subdomain.numprocs; i++)
        free(table[i].elms);
    free(table);
}

void _fmd_field_call_update_event_handler(fmd_t *md, int field_index, int turi_index)
{
    fmd_event_params_field_update_t params;

    params.field = field_index;
    params.turi = turi_index;

    md->EventHandler(md, FMD_EVENT_FIELD_UPDATE, (fmd_params_t *)&params);
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

static void gather_field_data_real(fmd_t *md, turi_t *t, field_t *f, fmd_array3s_t *out)
{
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* only "owners" participate in the collective communication below */
    if (owcomm->owned_tcells_num == 0) return;

    fmd_real_t ***lcd = (fmd_real_t ***)f->data.data;

    fmd_real_t *sendbuf = (fmd_real_t *)m_alloc(owcomm->owned_tcells_num * sizeof(*sendbuf));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        sendbuf[i] = ARRAY_ELEMENT(lcd, itc);
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
            ARRAY_ELEMENT((fmd_real_t ***)out->data, itc) = recvbuf[i];
        }

        free(recvbuf);
    }
}

static void gather_field_data_unsigned(fmd_t *md, turi_t *t, field_t *f, fmd_array3s_t *out)
{
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* only "owners" participate in the collective communication below */
    if (owcomm->owned_tcells_num == 0) return;

    unsigned ***lcd = (unsigned ***)f->data.data;

    unsigned *sendbuf = (unsigned *)m_alloc(owcomm->owned_tcells_num * sizeof(*sendbuf));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        sendbuf[i] = ARRAY_ELEMENT(lcd, itc);
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
            ARRAY_ELEMENT((unsigned ***)out->data, itc) = recvbuf[i];
        }

        free(recvbuf);
    }
}

static void gather_field_data_rtuple(fmd_t *md, turi_t *t, field_t *f, fmd_array3s_t *out)
{
    turi_ownerscomm_t *owcomm = &t->ownerscomm;

    /* only "owners" participate in the collective communication below */
    if (owcomm->owned_tcells_num == 0) return;

    fmd_rtuple_t ***lcd = (fmd_rtuple_t ***)f->data.data;

    fmd_rtuple_t *sendbuf = (fmd_rtuple_t *)m_alloc(owcomm->owned_tcells_num * sizeof(*sendbuf));

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        for (int d=0; d<DIM; d++)
            sendbuf[i][d] = ARRAY_ELEMENT(lcd, itc)[d];
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

            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT((fmd_rtuple_t ***)out->data, itc)[d] = recvbuf[i][d];
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

        if (pset[0] == md->Subdomain.myrank)
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

/* Finds in each direction the ranks of the upper and lower subdomains which are "owners" of
   turi-cells, and initialize rank_of_lower_owner and rank_of_upper_owner accordingly. */
static void init_rank_of_lower_upper_owner(fmd_t *md, turi_t *t)
{
    fmd_ituple_t is_tempo;

    for (int d=0; d<DIM; d++)                 /* set is_tempo[] to current subdomain */
        is_tempo[d] = md->Subdomain.is[d];

    for (int d=0; d<DIM; d++)
    {
        fmd_real_t pos;

        /* initialize rank_of_upper_owner[d] */

        pos = (t->itc_stop_global[d] % t->tdims_global[d]) * t->tcellh[d];
        is_tempo[d] = (int)impreal( _fmd_convert_pos_to_subd_coord_1D(md, pos, d) );
        t->rank_of_upper_owner[d] = INDEX_FLAT(is_tempo, md->ns);

        if ( (t->rank_of_upper_owner[d] < md->Subdomain.myrank && !md->PBC[d]) ||
             (t->rank_of_upper_owner[d] == md->Subdomain.myrank) )
            t->rank_of_upper_owner[d] = MPI_PROC_NULL;

        /* initialize rank_of_lower_owner[d] */

        int ifst = t->itc_start_owned[d] - t->itc_glob_to_loc[d]; /* convert to global index */
        pos = ((ifst - 1 + t->tdims_global[d]) % t->tdims_global[d]) * t->tcellh[d];
        is_tempo[d] = (int)impreal( _fmd_convert_pos_to_subd_coord_1D(md, pos, d) );
        t->rank_of_lower_owner[d] = INDEX_FLAT(is_tempo, md->ns);

        if ( (t->rank_of_lower_owner[d] > md->Subdomain.myrank && !md->PBC[d]) ||
             (t->rank_of_lower_owner[d] == md->Subdomain.myrank) )
            t->rank_of_lower_owner[d] = MPI_PROC_NULL;

        /* initialize has_upper_lower_procs[d] */

        t->has_upper_lower_owner_procs[d] = (t->rank_of_upper_owner[d] != MPI_PROC_NULL) ||
                                            (t->rank_of_lower_owner[d] != MPI_PROC_NULL);

        /* set is_tempo[] back again to current subdomain so as to prepare for the next iteration */

        is_tempo[d] = md->Subdomain.is[d];
    }
}

fmd_handle_t fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz, fmd_real_t starttime, fmd_real_t stoptime)
{
    if (md->Subdomain.grid == NULL) _fmd_subd_init(md);

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
        switch (cat)
        {
            case FMD_TURI_TTM_TYPE1:
                /* if ttm is 3D, all directions have margins, if it's 1D only one direction has margin */
                if (d == 2)
                    t->itc_start[d] = 1;
                else /* if (t->tdims_global[0] == 1 && t->tdims_global[1] == 1) is true, ttm is 1D */
                    t->itc_start[d] = (t->tdims_global[0] == 1 && t->tdims_global[1] == 1) ? 0 : 1;

                break;

            default:
                t->itc_start[d] = 0;
        }

        t->tcellh[d] = md->l[d] / t->tdims_global[d];

        fmd_real_t xlo = md->Subdomain.ic_global_firstcell[d] * md->cellh[d];
        fmd_real_t xlodiv = impreal(xlo / t->tcellh[d]);
        t->itc_start_global[d] = (int)xlodiv;

        /* if the fractional part of xlodiv is zero, then t->itc_start_owned[d] = t->itc_start[d] */
        t->itc_start_owned[d] = t->itc_start[d] + ((fmd_real_t)t->itc_start_global[d] == xlodiv ? 0 : 1);

        fmd_real_t xhi = xlo + md->Subdomain.cell_num_nonmarg[d] * md->cellh[d];
        t->itc_stop_global[d] = (int)ceil(impreal(xhi / t->tcellh[d]));

        t->tdims_local_nonmarg[d] = t->itc_stop_global[d] - t->itc_start_global[d];

        t->itc_stop[d] = t->itc_start[d] + t->tdims_local_nonmarg[d];

        t->tdims_local[d] = t->itc_start[d] + t->itc_stop[d];

        t->itc_glob_to_loc[d] = -t->itc_start_global[d] + t->itc_start[d];
    }

    t->tcell_volume = t->tcellh[0] * t->tcellh[1] * t->tcellh[2];

    t->cat = cat;

    t->ttm = NULL;
    t->fields = NULL;
    t->fields_num = 0;

    init_rank_of_lower_upper_owner(md, t);

    md->turies_num++;

    prepare_turi_for_communication(md, t);

    switch (cat)
    {
        case FMD_TURI_CUSTOM:
            break;

        case FMD_TURI_TTM_TYPE1:
            t->ttm = _fmd_ttm_construct(md, t);
            md->active_ttm_turi = t;
            break;

        default:
            assert(0);  /* TO-DO */
    }

    return ti;
}

static void set_field_data_el_size_and_type(field_t *f)
{
    if (f->cat == FMD_FIELD_MASS || f->cat == FMD_FIELD_TEMPERATURE ||
        f->cat == FMD_FIELD_NUMBER_DENSITY || f->cat == FMD_FIELD_TTM_TE ||
        f->cat == FMD_FIELD_TTM_XI)
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

inline static void prepare_buf1_for_reduce_real(fmd_real_t *buf1, turi_comm_t *tcm, field_t *f)
{
    for (int j=0; j < tcm->num_tcells; j++)
    {
        int *itc = tcm->itcs[j];

        buf1[j] = ARRAY_ELEMENT((fmd_real_t ***)f->data.data, itc);
    }
}

inline static void prepare_buf1_for_reduce_rtuple(fmd_rtuple_t *buf1, turi_comm_t *tcm, field_t *f)
{
    for (int j=0; j < tcm->num_tcells; j++)
    {
        int *itc = tcm->itcs[j];

        for (int d=0; d<DIM; d++)
            buf1[j][d] = ARRAY_ELEMENT((fmd_rtuple_t ***)f->data.data, itc)[d];
    }
}

inline static void prepare_buf1_for_reduce_unsigned(unsigned *buf1, turi_comm_t *tcm, field_t *f)
{
    for (int j=0; j < tcm->num_tcells; j++)
    {
        int *itc = tcm->itcs[j];

        buf1[j] = ARRAY_ELEMENT((unsigned ***)f->data.data, itc);
    }
}

static void perform_field_reduce_rtuple(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave)
{
    /* allocate memory for buf1 and buf2 */
    fmd_rtuple_t *buf1 = m_alloc(t->num_tcells_max * sizeof(*buf1));
    fmd_rtuple_t *buf2 = m_alloc(t->num_tcells_max * sizeof(*buf2));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buf1 (send-buffer) */

        prepare_buf1_for_reduce_rtuple(buf1, tcm, f);

        /* do communication */

        if (allhave)
            MPI_Allreduce(buf1, buf2, DIM*tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, tcm->comm);
        else
            MPI_Reduce(buf1, buf2, DIM*tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, 0, tcm->comm);

        /* copy from buf2 to data array */

        if (allhave || tcm->pset[0] == md->Subdomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                for (int d=0; d<DIM; d++)
                    ARRAY_ELEMENT((fmd_rtuple_t ***)f->data.data, itc)[d] = ((fmd_rtuple_t *)buf2)[j][d];
            }
        }
    }

    free(buf1);
    free(buf2);
}

static void perform_field_reduce_real(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave)
{
    /* allocate memory for buf1 and buf2 */
    fmd_real_t *buf1 = m_alloc(t->num_tcells_max * sizeof(*buf1));
    fmd_real_t *buf2 = m_alloc(t->num_tcells_max * sizeof(*buf2));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buf1 (send-buffer) */

        prepare_buf1_for_reduce_real(buf1, tcm, f);

        /* do communication */

        if (allhave)
            MPI_Allreduce(buf1, buf2, tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, tcm->comm);
        else
            MPI_Reduce(buf1, buf2, tcm->num_tcells, FMD_MPI_REAL, MPI_SUM, 0, tcm->comm);

        /* copy from buf2 to data array */

        if (allhave || tcm->pset[0] == md->Subdomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                ARRAY_ELEMENT((fmd_real_t ***)f->data.data, itc) = ((fmd_real_t *)buf2)[j];
            }
        }
    }

    free(buf1);
    free(buf2);
}

static void perform_field_reduce_unsigned(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave)
{
    /* allocate memory for buf1 and buf2 */
    unsigned *buf1 = m_alloc(t->num_tcells_max * sizeof(*buf1));
    unsigned *buf2 = m_alloc(t->num_tcells_max * sizeof(*buf2));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buf1 (send-buffer) */

        prepare_buf1_for_reduce_unsigned(buf1, tcm, f);

        /* do communication */

        if (allhave)
            MPI_Allreduce(buf1, buf2, tcm->num_tcells, MPI_UNSIGNED, MPI_SUM, tcm->comm);
        else
            MPI_Reduce(buf1, buf2, tcm->num_tcells, MPI_UNSIGNED, MPI_SUM, 0, tcm->comm);

        /* copy from buf2 to data array */

        if (allhave || tcm->pset[0] == md->Subdomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                ARRAY_ELEMENT((unsigned ***)f->data.data, itc) = ((unsigned *)buf2)[j];
            }
        }

    }

    free(buf1);
    free(buf2);
}

static void perform_field_bcast_real(fmd_t *md, field_t *f, turi_t *t)
{
    /* allocate memory for buffer */
    fmd_real_t *buf = m_alloc(t->num_tcells_max * sizeof(*buf));

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        /* prepare buffer on root process */

        if (tcm->pset[0] == md->Subdomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                buf[j] = ARRAY_ELEMENT((fmd_real_t ***)f->data.data, itc);
            }
        }

        /* do communication (broadcast buffer) */

        MPI_Bcast(buf, tcm->num_tcells, FMD_MPI_REAL, 0, tcm->comm);

        /* copy from buffer to data array on non-root processes */

        if (tcm->pset[0] != md->Subdomain.myrank)
        {
            for (int j=0; j < tcm->num_tcells; j++)
            {
                int *itc = tcm->itcs[j];

                ARRAY_ELEMENT((fmd_real_t ***)f->data.data, itc) = ((fmd_real_t *)buf)[j];
            }
        }
    }

    free(buf);
}

static void update_field_number(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                                int time_iteration, fmd_real_t time)
{
    unsigned ***num = (unsigned ***)f->data.data;

    /* clean data of number field (initialize with zeros) */
    _fmd_array_3d_unsigned_clean(num, t->tdims_local);

    fmd_ituple_t ic, itc;
    cell_t *cell;
    int i;
    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < cell->parts_num; i++)
        {
            /* find local index of turi-cell from particle's position */
            for (int d=0; d<DIM; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* count atoms */
            ARRAY_ELEMENT(num, itc)++;
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_reduce_unsigned(md, f, t, allhave);

    f->timestep = time_iteration; /* mark as updated */
    f->time = time;

    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_number_density(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                                        int time_iteration, fmd_real_t time)
{
    field_t *fdep = &t->fields[f->dependcs[0]];

    /* update the dependency field if not already updated */
    if (fdep->timestep != time_iteration)
        update_field_number(md, fdep, t, allhave || fdep->allhave_now, time_iteration, time);

    unsigned ***num = (unsigned ***)fdep->data.data;
    fmd_real_t ***nd = (fmd_real_t ***)f->data.data;

    if (allhave)
    {
        fmd_ituple_t itc;

        LOOP3D(itc, t->itc_start, t->itc_stop)
            ARRAY_ELEMENT(nd, itc) = ARRAY_ELEMENT(num, itc) / t->tcell_volume;
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            ARRAY_ELEMENT(nd, itc) = ARRAY_ELEMENT(num, itc) / t->tcell_volume;
        }
    }

    f->timestep = time_iteration; /* mark as updated */
    f->time = time;

    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_mass(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                              int time_iteration, fmd_real_t time)
{
    fmd_real_t ***mass = (fmd_real_t ***)f->data.data;

    /* clean data of mass field (initialize with zeros) */
    _fmd_array_3d_real_clean(mass, t->tdims_local);

    fmd_ituple_t ic, itc;
    cell_t *cell;
    int i;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < cell->parts_num; i++)
        {
            /* find index of turi-cell from particle's position */
            for (int d=0; d<DIM; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* update mass for turi-cell with index of itc */
            ARRAY_ELEMENT(mass, itc) += md->potsys.atomkinds[cell->atomkind[i]].mass;
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_reduce_real(md, f, t, allhave);

    f->timestep = time_iteration; /* mark as updated */
    f->time = time;

    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_vcm_only(fmd_t *md, field_t *fvcm, field_t *fmass, turi_t *t, fmd_bool_t allhave,
                                  int time_iteration, fmd_real_t time)
{
    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)fvcm->data.data;

    /* clean data of vcm field (initialize with zeros) */
    _fmd_array_3d_rtuple_clean(vcm, t->tdims_local);

    cell_t *cell;
    int i;
    fmd_ituple_t ic, itc;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < cell->parts_num; i++)
        {

            /* find index of turi-cell from particle's position */
            for (int d=0; d<DIM; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* calculate vcm (first phase) */

            fmd_real_t m = md->potsys.atomkinds[cell->atomkind[i]].mass;

            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT(vcm, itc)[d] += m * VEL(cell, i, d);
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_reduce_rtuple(md, fvcm, t, allhave);

    /* calculate vcm (last phase) */

    fmd_real_t ***mass = (fmd_real_t ***)fmass->data.data;

    if (allhave)
    {
        LOOP3D(itc, t->itc_start, t->itc_stop)
            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT(vcm, itc)[d] /= ARRAY_ELEMENT(mass, itc);
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT(vcm, itc)[d] /= ARRAY_ELEMENT(mass, itc);
        }
    }

    fvcm->timestep = time_iteration; /* mark as updated */
    fvcm->time = time;

    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, fvcm->field_index, t->turi_index);
}

static void update_field_vcm_and_mass(fmd_t *md, field_t *fvcm, field_t *fmass, turi_t *t, fmd_bool_t allhave,
                                      int time_iteration, fmd_real_t time)
{
    fmd_ituple_t ic, itc;

    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)fvcm->data.data;
    fmd_real_t ***mass = (fmd_real_t ***)fmass->data.data;

    /* clean data of mass and vcm fields (initialize with zeros) */
    LOOP3D(itc, t->itc_start, t->itc_stop)
    {
        ARRAY_ELEMENT(mass, itc) = 0.0;
        for (int d=0; d<DIM; d++)
            ARRAY_ELEMENT(vcm, itc)[d] = 0.0;
    }

    cell_t *cell;
    int i;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < cell->parts_num; i++)
        {

            /* find index of turi-cell from particle's position */
            for (int d=0; d<DIM; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* calculate vcm and mass (first phase) */

            fmd_real_t m = md->potsys.atomkinds[cell->atomkind[i]].mass;

            ARRAY_ELEMENT(mass, itc) += m;
            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT(vcm, itc)[d] += m * VEL(cell, i, d);
        }

    /* do communications */
    if (t->comms_num > 0)
    {
        perform_field_reduce_real(md, fmass, t, allhave || fmass->allhave_now);
        perform_field_reduce_rtuple(md, fvcm, t, allhave);
    }

    /* calculate vcm (last phase) */
    if (allhave)
    {
        LOOP3D(itc, t->itc_start, t->itc_stop)
            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT(vcm, itc)[d] /= ARRAY_ELEMENT(mass, itc);
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            for (int d=0; d<DIM; d++)
                ARRAY_ELEMENT(vcm, itc)[d] /= ARRAY_ELEMENT(mass, itc);
        }
    }

    fmass->timestep = time_iteration; /* mark as updated */
    fmass->time = time;
    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, fmass->field_index, t->turi_index);

    fvcm->timestep = time_iteration; /* mark as updated */
    fvcm->time = time;
    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, fvcm->field_index, t->turi_index);
}

static void update_field_vcm(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                             int time_iteration, fmd_real_t time)
{
    field_t *fmass = &t->fields[f->dependcs[0]];

    if (fmass->timestep != time_iteration)  /* if fmass is not already updated */
        update_field_vcm_and_mass(md, f, fmass, t, allhave, time_iteration, time);
    else
        update_field_vcm_only(md, f, fmass, t, allhave, time_iteration, time);
}

static void update_field_temperature_and_number(fmd_t *md, field_t *ftemp, field_t *fnum, field_t *fvcm,
                                                turi_t *t, fmd_bool_t allhave, int time_iteration,
                                                fmd_real_t time)
{
    fmd_real_t ***temp = (fmd_real_t ***)ftemp->data.data;
    unsigned ***num = (unsigned ***)fnum->data.data;
    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)fvcm->data.data;

    fmd_ituple_t itc;

    /* clean data of temperature and number fields (initialize with zeros) */
    LOOP3D(itc, t->itc_start, t->itc_stop)
    {
        ARRAY_ELEMENT(temp, itc) = 0.0;
        ARRAY_ELEMENT(num, itc) = 0;
    }

    cell_t *cell;
    int i;
    fmd_ituple_t ic;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < cell->parts_num; i++)
        {
            /* find index of turi-cell from particle's position */
            for (int d=0; d<DIM; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* calculate temperature (first phase) */

            fmd_rtuple_t vt;

            diffrt(vt, &VEL(cell, i, 0), ARRAY_ELEMENT(vcm, itc));

            fmd_real_t m = md->potsys.atomkinds[cell->atomkind[i]].mass;

            ARRAY_ELEMENT(temp, itc) += m * sqrrt(vt);

            /* count the particle in number field */
            ++ARRAY_ELEMENT(num, itc);
        }

    /* do communications */
    if (t->comms_num > 0)
    {
        perform_field_reduce_unsigned(md, fnum, t, allhave || fnum->allhave_now);
        perform_field_reduce_real(md, ftemp, t, allhave);
    }

    /* calculate temperature (last phase) */

    if (allhave)
    {
        LOOP3D(itc, t->itc_start, t->itc_stop)
            ARRAY_ELEMENT(temp, itc) /= (3.0 * K_BOLTZMANN * ARRAY_ELEMENT(num, itc));
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            ARRAY_ELEMENT(temp, itc) /= (3.0 * K_BOLTZMANN * ARRAY_ELEMENT(num, itc));
        }
    }

    fnum->timestep = time_iteration; /* mark as updated */
    fnum->time = time;
    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, fnum->field_index, t->turi_index);

    ftemp->timestep = time_iteration; /* mark as updated */
    ftemp->time = time;
    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, ftemp->field_index, t->turi_index);
}

static void update_field_temperature_only(fmd_t *md, field_t *ftemp, field_t *fnum, field_t *fvcm,
                                          turi_t *t, fmd_bool_t allhave, int time_iteration,
                                          fmd_real_t time)
{
    fmd_real_t ***temp = (fmd_real_t ***)ftemp->data.data;
    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)fvcm->data.data;

    /* clean data of temperature field (initialize with zeros) */
    _fmd_array_3d_real_clean(temp, t->tdims_local);

    cell_t *cell;
    int i;
    fmd_ituple_t ic, itc;

    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < cell->parts_num; i++)
        {
            /* find index of turi-cell from particle's position */
            for (int d=0; d<DIM; d++)
                itc[d] = (int)(POS(cell, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

            /* calculate temperature (first phase) */

            fmd_rtuple_t vt;

            diffrt(vt, &VEL(cell, i, 0), ARRAY_ELEMENT(vcm, itc));

            fmd_real_t m = md->potsys.atomkinds[cell->atomkind[i]].mass;

            ARRAY_ELEMENT(temp, itc) += m * sqrrt(vt);
        }

    /* do communications */
    if (t->comms_num > 0) perform_field_reduce_real(md, ftemp, t, allhave);

    /* calculate temperature (last phase) */

    unsigned ***num = (unsigned ***)fnum->data.data;

    if (allhave)
    {
        LOOP3D(itc, t->itc_start, t->itc_stop)
            ARRAY_ELEMENT(temp, itc) /= (3.0 * K_BOLTZMANN * ARRAY_ELEMENT(num, itc));
    }
    else
    {
        for (int i=0; i < t->ownerscomm.owned_tcells_num; i++)
        {
            int *itc = t->ownerscomm.owned_tcells[i];
            ARRAY_ELEMENT(temp, itc) /= (3.0 * K_BOLTZMANN * ARRAY_ELEMENT(num, itc));
        }
    }

    ftemp->timestep = time_iteration; /* mark as updated */
    ftemp->time = time;

    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, ftemp->field_index, t->turi_index);
}

static void update_field_temperature(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                                     int time_iteration, fmd_real_t time)
{
    field_t *fnum = &t->fields[f->dependcs[0]];
    field_t *fvcm = &t->fields[f->dependcs[1]];

    /* update the vcm field if not already updated */
    if (fvcm->timestep != time_iteration)
        update_field_vcm(md, fvcm, t, FMD_TRUE, time_iteration, time);

    if (fnum->timestep != time_iteration)  /* if fnum is not already updated */
        update_field_temperature_and_number(md, f, fnum, fvcm, t, allhave, time_iteration, time);
    else
        update_field_temperature_only(md, f, fnum, fvcm, t, allhave, time_iteration, time);
}

static void update_field_ttm_Te_and_xi(fmd_t *md, field_t *f, turi_t *t, int time_iteration, fmd_real_t time)
{
    field_t *ftemp = &t->fields[f->dependcs[0]];

    /* update the temperature field if not already updated */
    if (ftemp->timestep != time_iteration)
        update_field_temperature(md, ftemp, t, FMD_FALSE, time_iteration, time);

    ttm_t *ttm = t->ttm;

    if (t->ownerscomm.owned_tcells_num > 0) ttm->update_xe_te(md, t, ttm);

    field_t *fxi = &t->fields[ttm->ixi];
    fxi->time = time;                    /* mark as updated */
    fxi->timestep = time_iteration;

    field_t *fTe = &t->fields[ttm->iTe];
    fTe->time = time;                    /* mark as updated */
    fTe->timestep = time_iteration;

    if (t->comms_num > 0) perform_field_bcast_real(md, fxi, t);

    if (md->EventHandler != NULL) _fmd_field_call_update_event_handler(md, fTe->field_index, t->turi_index);
}

#define ADD_interval_TO_ARRAY(array, num, interval)                            \
    do                                                                         \
    {                                                                          \
        (array) = (fmd_real_t *)re_alloc(array, ((num)+1)*sizeof(fmd_real_t)); \
        (array)[(num)++] = (interval);                                         \
    } while(0)

static void field_intervals_add(field_t *f, fmd_real_t interval, fmd_bool_t allhave)
{
    unsigned j;

    /* see if the interval already exists in "intervals_allhave" array */
    for (j=0; j < f->intervals_allhave_num; j++)
        if (f->intervals_allhave[j] == interval) break;

    if (j == f->intervals_allhave_num) /* no, doesn't exist there */
    {
        /* see if the interval already exists in "intervals" array */
        for (j=0; j < f->intervals_num; j++)
            if (f->intervals[j] == interval) break;

        if (j == f->intervals_num) /* doesn't exist in "intervals" either */
        {
            if (allhave) /* add it to "intervals_allhave" */
                ADD_interval_TO_ARRAY(f->intervals_allhave, f->intervals_allhave_num, interval);
            else /* add it to "intervals" */
                ADD_interval_TO_ARRAY(f->intervals, f->intervals_num, interval);
        }
        else /* it exists in "intervals" array */
        {
            if (allhave)
            {
                /* remove it from "intervals" */
                f->intervals[j] = f->intervals[--f->intervals_num];
                f->intervals = re_alloc(f->intervals, f->intervals_num * sizeof(fmd_real_t));

                /* add it to "intervals_allhave" */
                ADD_interval_TO_ARRAY(f->intervals_allhave, f->intervals_allhave_num, interval);
            }
        }
    }
}

/* when a field is updated and allhave is equal to FMD_FALSE, in the current subdomain
   the field is only updated in those turi-cells which are "owned" by the current subdomain. */
int _fmd_field_add(turi_t *t, fmd_field_t cat, fmd_real_t interval, fmd_bool_t allhave)
{
    int i;
    field_t *f;
    unsigned dep1, dep2;

    /* add dependency fields */
    switch (cat)
    {
        case FMD_FIELD_TEMPERATURE:
            dep1 = _fmd_field_add(t, FMD_FIELD_NUMBER, interval, allhave);
            dep2 = _fmd_field_add(t, FMD_FIELD_VCM, interval, FMD_TRUE);
            break;

        case FMD_FIELD_VCM:
            dep1 = _fmd_field_add(t, FMD_FIELD_MASS, interval, allhave);
            break;

        case FMD_FIELD_NUMBER_DENSITY:
            dep1 = _fmd_field_add(t, FMD_FIELD_NUMBER, interval, allhave);
            break;

        case FMD_FIELD_TTM_TE:
            dep1 = _fmd_field_add(t, FMD_FIELD_TEMPERATURE, interval, FMD_FALSE);
            break;

        case FMD_FIELD_TTM_XI:
            dep1 = _fmd_field_add(t, FMD_FIELD_TEMPERATURE, interval, FMD_FALSE);
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
        f->intervals_allhave = NULL;
        f->intervals_allhave_num = 0;
        f->intervals = NULL;
        f->intervals_num = 0;
        t->fields_num++;

        if (cat == FMD_FIELD_TEMPERATURE)
        {
            f->dependcs_num = 2;
            f->dependcs = (unsigned *)m_alloc(f->dependcs_num * sizeof(unsigned));
            f->dependcs[0] = dep1;
            f->dependcs[1] = dep2;
        }
        else if (cat == FMD_FIELD_VCM || cat == FMD_FIELD_NUMBER_DENSITY || cat == FMD_FIELD_TTM_TE ||
                 cat == FMD_FIELD_TTM_XI)
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

    /* add the interval if doesn't already exist in "intervals" or "intervals_allhave" arrays */
    field_intervals_add(f, interval, allhave);

    return i;
}

fmd_handle_t fmd_field_add(fmd_t *md, fmd_handle_t turi, fmd_field_t cat, fmd_real_t interval)
{
    return _fmd_field_add(&md->turies[turi], cat, interval, FMD_FALSE);
}

static void update_field_anykind(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                                 int time_iteration, fmd_real_t time)
{
    switch (f->cat)
    {
        case FMD_FIELD_NUMBER:
            update_field_number(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_NUMBER_DENSITY:
            update_field_number_density(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_MASS:
            update_field_mass(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_VCM:
            update_field_vcm(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_TEMPERATURE:
            update_field_temperature(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_TTM_TE:
        case FMD_FIELD_TTM_XI:
            update_field_ttm_Te_and_xi(md, f, t, time_iteration, time);
            break;
    }
}

static void update_field_ForceIndependent(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                                          int time_iteration, fmd_real_t time)
{
    switch (f->cat)
    {
        case FMD_FIELD_NUMBER:
            update_field_number(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_NUMBER_DENSITY:
            update_field_number_density(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_MASS:
            update_field_mass(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_VCM:
            update_field_vcm(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_TEMPERATURE:
            update_field_temperature(md, f, t, allhave, time_iteration, time);
            break;

        case FMD_FIELD_TTM_TE:
        case FMD_FIELD_TTM_XI:
            update_field_ttm_Te_and_xi(md, f, t, time_iteration, time);
            break;
    }
}

static void update_field(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allhave,
                         int time_iteration, fmd_real_t time,
                         fmd_bool_t Xupd, fmd_bool_t Vupd, fmd_bool_t Fupd)
{
    if (Xupd && Vupd && Fupd)
        update_field_anykind(md, f, t, allhave, time_iteration, time);
    else
        if (Xupd && Vupd)
            update_field_ForceIndependent(md, f, t, allhave, time_iteration, time);
        else
            assert(0); /* TO-DO - this kind is not handled yet */
}

/* this subroutine updates the variable allhave_now in all fields */
static void allhave_now_update(fmd_t *md, turi_t *t, fmd_bool_t IsTTM)
{
    for (int fi = 0; fi < t->fields_num; fi++)
    {
        field_t *f = &t->fields[fi];

        f->allhave_now = FMD_FALSE;

        /* fields on a TTM turi must be updated every time step */
        if (IsTTM && f->intervals_allhave_num > 0)
        {
            f->allhave_now = FMD_TRUE;
            continue;
        }

        for (int j = 0; j < f->intervals_allhave_num; j++)
            if (_fmd_timer_is_its_time(md->time, md->timestep/2.0, t->starttime, f->intervals_allhave[j]))
            {
                f->allhave_now = FMD_TRUE;
                break;
            }
    }
}

void _fmd_turies_update(fmd_t *md, int time_iteration, fmd_real_t time,
                        fmd_bool_t Xupd, fmd_bool_t Vupd, fmd_bool_t Fupd)
{
    for (int ti=0; ti < md->turies_num; ti++)
    {
        turi_t *t = &md->turies[ti];

        if (_is_time_within_turi_start_stop_times(md, t))
        {
            fmd_bool_t IsTTM = (t->cat == FMD_TURI_TTM_TYPE1);

            allhave_now_update(md, t, IsTTM);

            /* start from fields with higher indexes */
            for (int fi = t->fields_num-1; fi >= 0; fi--)
            {
                field_t *f = &t->fields[fi];

                if (f->timestep == time_iteration) continue; /* see if the field is already updated */

                if (f->allhave_now)
                    update_field(md, f, t, FMD_TRUE, time_iteration, time, Xupd, Vupd, Fupd);
                else
                {
                    if (IsTTM) /* fields on a TTM turi must be updated every time step */
                        update_field(md, f, t, FMD_FALSE, time_iteration, time, Xupd, Vupd, Fupd);
                    else
                    {
                        for (int j=0; j < f->intervals_num; j++)
                            if (_fmd_timer_is_its_time(md->time, md->timestep/2.0, t->starttime, f->intervals[j]))
                            {
                                update_field(md, f, t, FMD_FALSE, time_iteration, time, Xupd, Vupd, Fupd);
                                break;
                            }
                    }
                }
            }
        }
    }
}

static void convert_inmass_to_outmass(fmd_array3s_t *ar)
{
    fmd_utriple_t iv;

    LOOP3D(iv, _fmd_ThreeZeros_int, ar->dims)
        ARRAY_ELEMENT((fmd_real_t ***)(ar->data), iv) *= MD_MASS_UNIT;
}

fmd_array3s_t *fmd_field_getArray(fmd_t *md, fmd_handle_t turi, fmd_handle_t field,
                                  fmd_array3_t *array, fmd_utriple_t dims)
{
    turi_t *t = &md->turies[turi];
    field_t *f = &t->fields[field];

    fmd_array3s_t *p;

    if (md->Is_MD_comm_root) p = m_alloc(sizeof(*p));

    switch(f->cat)
    {
        case FMD_FIELD_VCM:
            gather_field_data_rtuple(md, t, f, p);
            break;

        case FMD_FIELD_NUMBER_DENSITY:
        case FMD_FIELD_TEMPERATURE:
        case FMD_FIELD_TTM_TE:
            gather_field_data_real(md, t, f, p);
            break;

        case FMD_FIELD_MASS:
            gather_field_data_real(md, t, f, p);
            if (md->Is_MD_comm_root) convert_inmass_to_outmass(p);
            break;

        case FMD_FIELD_NUMBER:
            gather_field_data_unsigned(md, t, f, p);
            break;

        default:
            assert(0);
    }

    if (md->Is_MD_comm_root)
    {
        *array = p->data;
        dims[0] = p->dims[0];
        dims[1] = p->dims[1];
        dims[2] = p->dims[2];
    }
    else
        p = NULL;

    return p;
}

void fmd_field_save_as_hdf5(fmd_t *md, fmd_handle_t turi, fmd_handle_t field, fmd_string_t path)
{
    turi_t *t = &md->turies[turi];
    field_t *f = &t->fields[field];

    fmd_array3s_t data;

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

        case FMD_FIELD_TEMPERATURE:
            gather_field_data_real(md, t, f, &data);
            if (md->Is_MD_comm_root)
            {
                _fmd_h5_save_scalar_field_float(md, "temperature", t, path, &data);
                _fmd_array_3d_free(&data);
            }
            break;

        case FMD_FIELD_TTM_TE:
            gather_field_data_real(md, t, f, &data);
            if (md->Is_MD_comm_root)
            {
                _fmd_h5_save_scalar_field_float(md, "ttm-Te", t, path, &data);
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
    free(f->intervals_allhave);
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

    if (t->ttm != NULL) _fmd_ttm_destruct(&t->ttm);

    turi_ownerscomm_free(md, &t->ownerscomm);
}

void fmd_turi_free(fmd_t *md)
{
    for (unsigned u=0; u < md->turies_num; u++)
        turi_free(md, &md->turies[u]);
    free(md->turies);

    md->turies = NULL;
    md->turies_num = 0;
    md->active_ttm_turi = NULL;
}

/* negative returned value means the field doesn't exist */
fmd_handle_t fmd_field_find(fmd_t *md, fmd_handle_t turi, fmd_field_t cat)
{
    turi_t *t = &md->turies[turi];

    for (unsigned u=0; u < t->fields_num; u++)
        if (t->fields[u].cat == cat) return u;

    return -1;
}
