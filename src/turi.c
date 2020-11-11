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
        _fmd_array_3d_create(t->tdims_global[0], t->tdims_global[1], t->tdims_global[2],
                             sizeof(unsigned), out);
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
        is_start[d] = (int)slo[d];
        is_stop[d] = (int)ceil(shi[d]);
        np *= is_stop[d] - is_start[d];
    }

    *pset = (int *)malloc(np * sizeof(int));
    /* TO-DO: handle memory error */
    assert(*pset != NULL);

    fmd_ituple_t is;
    int i=0;

    ITERATE(is, is_start, is_stop)
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

    ITERATE(itc, t->tcell_start, t->tcell_stop)
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

    t->tcells_global_num = dimx * dimy * dimz;

    for (int d=0; d<3; d++)
    {
        t->tcellh[d] = md->l[d] / t->tdims_global[d];

        fmd_real_t xlo = md->SubDomain.ic_global_firstcell[d] * md->cellh[d];
        t->tcell_start[d] = (int)(xlo / t->tcellh[d]);
        fmd_real_t xhi = xlo + md->SubDomain.cell_num_nonmarg[d] * md->cellh[d];
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
        default:
            assert(0);  /* TO-DO */
    }

    md->turies_num++;

    prepare_turi_for_communication(md, t);

    /* only for test */
    /* for (int i=0; i < md->ns[0]; i++)
        for (int j=0; j < md->ns[1]; j++)
            for (int k=0; k < md->ns[2]; k++)
            {
                MPI_Barrier(md->MD_comm);
                if (md->SubDomain.is[0] == i && md->SubDomain.is[1] == j && md->SubDomain.is[2] == k)
                    printf("subdomain(%d, %d, %d)\n", i, j, k);
            }
    */

    return ti;
}

static void set_field_data_el_size_and_type(field_t *f)
{
    if (f->cat == FMD_FIELD_MASS || f->cat == FMD_FIELD_TEMPERATURE)
    {
        f->data_el_size = sizeof(fmd_real_t);
        f->data_type = FIELD_DATA_REAL;
    }
    else if (f->cat == FMD_FIELD_VCM)
    {
        f->data_el_size = sizeof(fmd_rtuple_t);
        f->data_type = FIELD_DATA_RTUPLE;
    }
    else if (f->cat == FMD_FIELD_NUMBER)
    {
        f->data_el_size = sizeof(unsigned);
        f->data_type = FIELD_DATA_UNSIGNED;
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

static void update_field_number(fmd_t *md, field_t *f, turi_t *t)
{
    fmd_ituple_t ic, itc;

    unsigned ***num = (unsigned ***)f->data.data;

    /* clean data of number field (initialize with zeros) */
    _fmd_array_3d_unsigned_clean(num, t->tdims[0], t->tdims[1], t->tdims[2]);

    cell_t *cell;
    int i;
    /* iterate over all particles in current subdomain */
    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {

            /* find index of turi-cell from particle's position */
            for (int d=0; d<3; d++)
                itc[d] = (int)(cell->parts[i].core.x[d] / t->tcellh[d]) - t->tcell_start[d];

            /* count atoms */
            num[itc[0]][itc[1]][itc[2]]++;
        }

    /* do communications */
    perform_field_comm_unsigned(md, f, t, FMD_FALSE);
}

static void update_field_vcm(fmd_t *md, field_t *f, turi_t *t)
{
    fmd_ituple_t ic, itc;

    fmd_rtuple_t ***vcm = (fmd_rtuple_t ***)f->data.data;
    fmd_real_t ***mass = (fmd_real_t ***)t->fields[f->dependcs[0]].data.data;

    /* clean data of mass and vcm fields (initialize with zeros) */
    ITERATE(itc, _fmd_ThreeZeros, t->tdims)
    {
        mass[itc[0]][itc[1]][itc[2]] = 0.0;
        for (int d=0; d<3; d++)
            vcm[itc[0]][itc[1]][itc[2]][d] = 0.0;
    }

    cell_t *cell;
    int i;

    /* iterate over all particles in current subdomain */
    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
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
    ITERATE(itc, _fmd_ThreeZeros, t->tdims)
        for (int d=0; d<3; d++)
            vcm[itc[0]][itc[1]][itc[2]][d] /= mass[itc[0]][itc[1]][itc[2]];
}

/*
void fmd_field_update_test(fmd_t *md)
{
    turi_t *t = &md->turies[0];
    field_t *f = &t->fields[0];

    update_field_number(md, f, t);

    fmd_array3D_t num;

    gather_field_data_unsigned(md, t, f, &num);

    if (md->Is_MD_comm_root)
    {
        for (int i=0; i < num.dim1; i++)
            for (int j=0; j < num.dim2; j++)
                for (int k=0; k < num.dim3; k++)
                    printf("%d ", ((unsigned ***)(num.data))[i][j][k]);
        printf("\n");

        _fmd_array_3d_free(&num);
    }
}*/

unsigned fmd_field_add(fmd_t *md, unsigned turi, fmd_field_t cat, fmd_real_t interval)
{
    unsigned i;
    field_t *f;
    unsigned dep1;
    turi_t *t = &md->turies[turi];

    /* add dependency fields */
    switch (cat)
    {
        case FMD_FIELD_TEMPERATURE:
            dep1 = fmd_field_add(md, turi, FMD_FIELD_VCM, interval);
            break;

        case FMD_FIELD_VCM:
            dep1 = fmd_field_add(md, turi, FMD_FIELD_MASS, interval);
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

        f = &t->fields[t->fields_num++];
        f->cat = cat;
        f->timestep = -1;
        f->intervals = NULL;
        f->intervals_num = 0;

        if (cat == FMD_FIELD_TEMPERATURE || cat == FMD_FIELD_VCM)
        {
            f->dependcs_num = 1;
            f->dependcs = (unsigned *)malloc(sizeof(unsigned));
            f->dependcs[0] = dep1;
        }
        else
            f->dependcs_num = 0;

        set_field_data_el_size_and_type(f);

        /* allocate space for field data */
        _fmd_array_3d_create(t->tdims[0], t->tdims[1], t->tdims[2], f->data_el_size, &f->data);
        assert(f->data.data != NULL);
    }
    else
        f = &t->fields[i];

    /* add the interval if doesn't already exist in the intervals array */
    unsigned j;
    for (j=0; j < f->intervals_num; j++)
        if (f->intervals[j] == interval) break;

    if (j == f->intervals_num)
    {
        f->intervals = (fmd_real_t *)realloc(f->intervals, (f->intervals_num+1) * sizeof(fmd_real_t));
        f->intervals[f->intervals_num++] = interval;
    }

    return i;
}

/*void fmd_field_report(fmd_t *md, unsigned turi)
{
    turi_t *t = &md->turies[turi];
    printf("number of added fields = %d\n\n", t->fields_num);
    for (int i=0; i < t->fields_num; i++)
    {
        printf("FIELD %d:\n", i);
        field_t *f = &t->fields[i];
        printf("cat = %d\n", f->cat);
        printf("timestep = %d\n", f->timestep);
        printf("number of dependency fields = %d\n", f->dependcs_num);
        printf("number of intervals = %d\n\n", f->intervals_num);
    }
}*/
