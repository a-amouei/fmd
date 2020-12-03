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

static void call_field_update_event_handler(fmd_t *md, int field_index, int turi_index)
{
    fmd_event_params_field_update_t params;

    params.field = field_index;
    params.turi = turi_index;

    md->EventHandler(md, FMD_EVENT_FIELD_UPDATE, &params);
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

    fmd_real_t *sendbuf = (fmd_real_t *)malloc(owcomm->owned_tcells_num * sizeof(*sendbuf));
    assert(sendbuf != NULL);
    /* TO-DO: handle memory error */

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        sendbuf[i] = lcd[itc[0]][itc[1]][itc[2]];
    }

    fmd_real_t *recvbuf;

    if (md->Is_MD_comm_root)
    {
        recvbuf = (fmd_real_t *)malloc(t->tcells_global_num * sizeof(*recvbuf));
        assert(recvbuf != NULL);
        /* TO-DO: handle memory error */
    }

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

    unsigned *sendbuf = (unsigned *)malloc(owcomm->owned_tcells_num * sizeof(*sendbuf));
    assert(sendbuf != NULL);
    /* TO-DO: handle memory error */

    for (int i=0; i < owcomm->owned_tcells_num; i++)
    {
        int *itc = owcomm->owned_tcells[i];
        sendbuf[i] = lcd[itc[0]][itc[1]][itc[2]];
    }

    unsigned *recvbuf;

    if (md->Is_MD_comm_root)
    {
        recvbuf = (unsigned *)malloc(t->tcells_global_num * sizeof(*recvbuf));
        assert(recvbuf != NULL);
        /* TO-DO: handle memory error */
    }

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

    *pset = (int *)malloc(np * sizeof(int));
    /* TO-DO: handle memory error */
    assert(*pset != NULL);

    fmd_ituple_t is;
    int i=0;

    LOOP3D(is, is_start, is_stop)
        (*pset)[i++] = INDEX(is, md->ns);

    /* if root process of the MD communicator exists in pset, let it
    occupy the first array element, so that it becomes an "owner". */
    for (int i=1; i < np; i++)
        if ((*pset)[i] == ROOTPROCESS(md->SubDomain.numprocs))
        {
            int process = (*pset)[0];
            (*pset)[0] = (*pset)[i];
            (*pset)[i] = process;
        }

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
    {
        owcomm->recvcounts = (int *)malloc(owcomm->commsize * sizeof(int));
        assert(owcomm->recvcounts != NULL);
        /* TO-DO: handle memory error */
    }

    MPI_Gather(&owcomm->owned_tcells_num, 1, MPI_INT, owcomm->recvcounts, 1, MPI_INT, 0, owcomm->comm);

    if (md->Is_MD_comm_root)
    {
        owcomm->global_indexes = (fmd_ituple_t *)malloc(t->tcells_global_num * sizeof(fmd_ituple_t));
        assert(owcomm->global_indexes != NULL);
        /* TO-DO: handle memory error */

        owcomm->displs = (int *)malloc(owcomm->commsize * sizeof(int));
        assert(owcomm->displs != NULL);
        /* TO-DO: handle memory error */

        owcomm->displs[0] = 0;
        for (int i=1; i < owcomm->commsize; i++)
            owcomm->displs[i] = owcomm->displs[i-1] + owcomm->recvcounts[i-1];
    }

    fmd_ituple_t *lg = (fmd_ituple_t *)malloc(owcomm->owned_tcells_num * sizeof(fmd_ituple_t));
    assert(lg != NULL);
    /* TO-DO: handle memory error */

    for (int i=0; i < owcomm->owned_tcells_num; i++)
        for (int d=0; d<3; d++)
            lg[i][d] = owcomm->owned_tcells[i][d] + t->tcell_start[d]; /* convert local indexes to global indexes */

    MPI_Gatherv(lg, owcomm->owned_tcells_num, md->mpi_types.mpi_ituple, owcomm->global_indexes,
                owcomm->recvcounts, owcomm->displs, md->mpi_types.mpi_ituple, 0, owcomm->comm);

    free(lg);
}

static void prepare_turi_for_communication(fmd_t *md, turi_t *t)
{
    t->comms_num = 0;
    t->comms = NULL;
    t->ownerscomm.owned_tcells = NULL;
    t->ownerscomm.owned_tcells_num = 0;

    MPI_Group mdgroup, newgroup;

    MPI_Comm_group(md->MD_comm, &mdgroup);

    fmd_ituple_t itc;

    LOOP3D(itc, t->tcell_start, t->tcell_stop)
    {
        turi_comm_t *tcomm;
        int *pset, np;

        np = identify_tcell_processes_set(md, t->tcellh, itc, &pset);

        /* do owned_tcells and owned_tcells_num need an update? */
        if (pset[0] == md->SubDomain.myrank)
        {
            t->ownerscomm.owned_tcells = (fmd_ituple_t *)realloc(t->ownerscomm.owned_tcells,
                                          (t->ownerscomm.owned_tcells_num+1) * sizeof(*t->ownerscomm.owned_tcells));
            /* TO-DO: handle memory error */
            assert(t->ownerscomm.owned_tcells != NULL);

            for (int d=0; d<3; d++)
                t->ownerscomm.owned_tcells[t->ownerscomm.owned_tcells_num][d] = itc[d] - t->tcell_start[d];

            t->ownerscomm.owned_tcells_num++;
        }

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
                MPI_Group_incl(mdgroup, np, pset, &newgroup);
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

        tcomm->itcs = (fmd_ituple_t *)realloc(tcomm->itcs, (tcomm->num_tcells+1) * sizeof(fmd_ituple_t));
        /* TO-DO: handle memory error */
        assert(tcomm->itcs != NULL);

        for (int d=0; d<3; d++)
            tcomm->itcs[tcomm->num_tcells][d] = itc[d] - t->tcell_start[d];

        tcomm->num_tcells++;
    }

    MPI_Group_free(&mdgroup);

    prepare_ownerscomm(md, t);

    /* find num_tcells_max */
    t->num_tcells_max = 0;
    for (int i=0; i < t->comms_num; i++)
        if (t->comms[i].num_tcells > t->num_tcells_max) t->num_tcells_max = t->comms[i].num_tcells;
}

fmd_handle_t fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz, fmd_real_t starttime, fmd_real_t stoptime)
{
    if (md->SubDomain.grid == NULL) fmd_subd_init(md);

    int ti = md->turies_num;

    md->turies = (turi_t *)realloc(md->turies, (ti+1) * sizeof(turi_t));
    /* TO-DO: handle memory error */
    assert(md->turies != NULL);

    turi_t *t = &md->turies[ti];

    t->turi_index = ti;

    t->tdims_global[0] = dimx;
    t->tdims_global[1] = dimy;
    t->tdims_global[2] = dimz;

    t->starttime = starttime;
    t->stoptime = stoptime;

    t->tcells_global_num = dimx * dimy * dimz;

    for (int d=0; d<3; d++)
    {
        t->tcellh[d] = md->l[d] / t->tdims_global[d];

        fmd_real_t xlo = md->SubDomain.ic_global_firstcell[d] * md->cellh[d];
        t->tcell_start[d] = (int)impreal(xlo / t->tcellh[d]);
        fmd_real_t xhi = xlo + md->SubDomain.cell_num_nonmarg[d] * md->cellh[d];
        t->tcell_stop[d] = (int)ceil(impreal(xhi / t->tcellh[d]));

        t->tdims[d] = t->tcell_stop[d] - t->tcell_start[d];
    }

    t->tcell_volume = t->tcellh[0] * t->tcellh[1] * t->tcellh[2];

    t->cat = cat;
    switch (cat)
    {
        case FMD_TURI_CUSTOM:
            t->fields = NULL;
            t->fields_num = 0;
            break;
        default:
            assert(0);  /* TO-DO */
    }

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

static void perform_field_comm_unsigned(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    /* allocate memory for buf1 and buf2 */
    unsigned *buf1 = malloc(t->num_tcells_max * sizeof(*buf1));
    assert(buf1 != NULL);
    unsigned *buf2 = malloc(t->num_tcells_max * sizeof(*buf2));
    assert(buf2 != NULL);
    /* TO-DO: handle memory error */

    /* perform communication for every communicator */

    for (int i=0; i < t->comms_num; i++)
    {
        turi_comm_t *tcm = &t->comms[i];

        if (tcm->commsize > 1)
        {
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
    }

    free(buf1);
    free(buf2);
}

static void update_field_number(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    fmd_ituple_t ic, itc;

    unsigned ***num = (unsigned ***)f->data.data;

    /* clean data of number field (initialize with zeros) */
    _fmd_array_3d_unsigned_clean(num, t->tdims);

    cell_t *cell;
    int i;
    /* iterate over all particles in current subdomain */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {
            /* find index of turi-cell from particle's position */
            for (int d=0; d<3; d++)
                itc[d] = (int)(cell->parts[i].core.x[d] / t->tcellh[d]) - t->tcell_start[d];

            /* count atoms */
            num[itc[0]][itc[1]][itc[2]]++;
        }

    /* do communications */
    perform_field_comm_unsigned(md, f, t, allreduce);

    f->timestep = md->time_iteration; /* mark as updated */

    if (md->EventHandler != NULL) call_field_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_number_density(fmd_t *md, field_t *f, turi_t *t, fmd_bool_t allreduce)
{
    field_t *fdep = &t->fields[f->dependcs[0]];
    /* update the dependency field if not already updated */
    if (fdep->timestep != md->time_iteration) update_field_number(md, fdep, t, allreduce);

    unsigned ***num = (unsigned ***)fdep->data.data;
    fmd_real_t ***nd = (fmd_real_t ***)f->data.data;

    fmd_ituple_t itc;

    LOOP3D(itc, _fmd_ThreeZeros_int, t->tdims)
        nd[itc[0]][itc[1]][itc[2]] = num[itc[0]][itc[1]][itc[2]] / t->tcell_volume;

    f->timestep = md->time_iteration; /* mark as updated */

    if (md->EventHandler != NULL) call_field_update_event_handler(md, f->field_index, t->turi_index);
}

static void update_field_vcm(fmd_t *md, field_t *f, turi_t *t)
{
    fmd_ituple_t ic, itc;

    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)f->data.data;
    fmd_real_t ***mass = (fmd_real_t ***)t->fields[f->dependcs[0]].data.data;

    /* clean data of mass and vcm fields (initialize with zeros) */
    LOOP3D(itc, _fmd_ThreeZeros_int, t->tdims)
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
                itc[d] = (int)(cell->parts[i].core.x[d] / t->tcellh[d]) - t->tcell_start[d];

            /* calculate vcm and mass (first phase) */

            fmd_real_t m = md->potsys.atomkinds[cell->parts[i].core.atomkind].mass;

            mass[itc[0]][itc[1]][itc[2]] += m;
            for (int d=0; d<3; d++)
                vcm[itc[0]][itc[1]][itc[2]][d] += m * cell->parts[i].core.v[d];
        }

    /* TO-DO: do communications */

    /* calculate vcm and mass (last phase) */
    LOOP3D(itc, _fmd_ThreeZeros_int, t->tdims)
        for (int d=0; d<3; d++)
            vcm[itc[0]][itc[1]][itc[2]][d] /= mass[itc[0]][itc[1]][itc[2]];
}

#define ADD_interval_TO_ARRAY(array, num, interval)                            \
    do                                                                         \
    {                                                                          \
        (array) = (fmd_real_t *)realloc(array, ((num)+1)*sizeof(fmd_real_t));  \
        assert((array) != NULL);      /* TO-DO */                              \
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
            dep1 = field_add(md, turi, FMD_FIELD_MASS, interval, FMD_TRUE);
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

        t->fields = (field_t *)realloc(t->fields, (t->fields_num+1) * sizeof(field_t));
        /* TO-DO: handle memory error */
        assert(t->fields != NULL);

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
            f->dependcs = (unsigned *)malloc(sizeof(unsigned));
            f->dependcs[0] = dep1;
        }
        else
            f->dependcs_num = 0;

        set_field_data_el_size_and_type(f);

        /* allocate space for field data */
        _fmd_array_3d_create(t->tdims, f->data_el_size, f->datatype, &f->data);
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
    }
}

void _fmd_turies_update(fmd_t *md)
{
    for (int ti=0; ti < md->turies_num; ti++)
    {
        turi_t *t = &md->turies[ti];

        if (md->MD_time >= t->starttime && !(md->MD_time > t->stoptime && t->stoptime >= t->starttime))
        {
            for (int fi=0; fi < t->fields_num; fi++)
            {
                field_t *f = &t->fields[fi];

                if (f->timestep == md->time_iteration) continue; /* see if it is already updated */

                int j;

                for (j=0; j < f->intervals_allreduce_num; j++)
                {
                    if (_fmd_timer_is_its_time(md->MD_time, md->delta_t/2.0, t->starttime, f->intervals_allreduce[j]))
                    {
                        update_field(md, f, t, FMD_TRUE);

                        break;
                    }
                }

                if (j == f->intervals_allreduce_num)
                {
                    for (j=0; j < f->intervals_num; j++)
                    {
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
}

void fmd_field_save_as_hdf5(fmd_t *md, fmd_handle_t turi, fmd_handle_t field, fmd_string_t path)
{
    turi_t *t = &md->turies[turi];
    field_t *f = &t->fields[field];

    fmd_array3D_t data;

    switch(f->cat)
    {
        case FMD_FIELD_NUMBER_DENSITY:
            gather_field_data_real(md, t, f, &data);
            if (md->Is_MD_comm_root)
            {
                _fmd_h5_save_scalar_float(md, "number-density", t, path, &data);
                _fmd_array_3d_free(&data);
            }
            break;

        default:
            assert(0);
    }
}
