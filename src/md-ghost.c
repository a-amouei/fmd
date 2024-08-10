/*
  md_ghost.c: This file is part of Free Molecular Dynamics

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

/* functions for working with ghost cells of MD subdomains, including
   communication of data from or to them */

#include <string.h>
#include "fmd-private.h"
#include "misc.h"
#include "md-ghost.h"
#include "types.h"
#include "general.h"
#include "cell.h"

/* *size is the size of the data in bytes, not the size of the buffer. The buffer can be larger. */
typedef void (*packer_t)(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                         int *size, fmd_pointer_t *out, bool nodest);
/* insize is the size of the buffer; dir is the direction of transfer (+1 or -1) */
typedef void (*unpacker_t)(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                           int insize, fmd_pointer_t in, int dim, int dir);

/* dir is the direction of transfer (+1 or -1) */
typedef void (*ccopier_t)(fmd_t *md, cell_t *src, cell_t *dest, int dim, int dir);

static unsigned set_partsnum_to_zero(fmd_t *md, fmd_ituple_t ic_from, fmd_ituple_t ic_to)
{
    fmd_ituple_t ic;
    unsigned count = 0;

    LOOP3D(ic, ic_from, ic_to)
    {
        cell_t *c = ARRAY_ELEMENT(md->subd.gridp, ic);

        count += c->parts_num;
        c->parts_num = 0;
    }

    return count;
}

void _fmd_ghostparticles_clean(fmd_t *md)
{
    for (int i = md->subd.nc; i < md->subd.ncm; i++)
        (md->subd.grid + i)->parts_num = 0;
}

static void ccopy_for_ghostinit(fmd_t *md, cell_t *src, cell_t *dest, int dim, int dir)
{
    dest->parts_num = src->parts_num;

    if (dest->capacity < dest->parts_num || dest->parts_num + md->cell_inc_dbl < dest->capacity)
        _fmd_cell_resize(md, dest);

    memcpy(dest->x, src->x, src->parts_num * sizeof(fmd_rtuple_t));
    memcpy(dest->GroupID, src->GroupID, src->parts_num * sizeof(int));
    memcpy(dest->atomkind, src->atomkind, src->parts_num * sizeof(unsigned));

    fmd_real_t value;

    if (dir > 0)
        value = -md->l[dim];
    else
        value = md->l[dim];

    for (int i=0; i < dest->parts_num; i++)
        POS(dest, i, dim) += value;
}

static void ccopy_for_partialforces(fmd_t *md, cell_t *src, cell_t *dest, int dim, int dir)
{
    for (unsigned i=0; i < src->parts_num; i++)
        for (int d=0; d < DIM; d++)
            FRC(dest, i, d) += FRC(src, i, d);
}

static void ccopy_for_Femb(fmd_t *md, cell_t *src, cell_t *dest, int dim, int dir)
{
    memcpy(dest->vaream, src->vaream, src->parts_num * sizeof(fmd_real_t));
}

static void ccopy_for_migrate(fmd_t *md, cell_t *src, cell_t *dest, int dim, int dir)
{
    int oldnum = dest->parts_num;

    dest->parts_num += src->parts_num;

    if (dest->capacity < dest->parts_num) _fmd_cell_resize(md, dest);

    memcpy(dest->x + DIM*oldnum, src->x, src->parts_num * sizeof(fmd_rtuple_t));

    if (dest->v != NULL)
        memcpy(dest->v + DIM*oldnum, src->v, src->parts_num * sizeof(fmd_rtuple_t));

    if (dest->GroupID != NULL)
        memcpy(dest->GroupID + oldnum, src->GroupID, src->parts_num * sizeof(int));

    if (dest->AtomID != NULL)
        memcpy(dest->AtomID + oldnum, src->AtomID, src->parts_num * sizeof(unsigned));

    memcpy(dest->atomkind + oldnum, src->atomkind, src->parts_num * sizeof(unsigned));

    fmd_real_t value;

    if (dir > 0)
        value = -md->l[dim];
    else
        value = md->l[dim];

    for (int i=oldnum; i < dest->parts_num; i++)
        POS(dest, i, dim) += value;

    src->parts_num = 0;
}

/* called instead of transfer() when ns==1 */
static void copy_local(fmd_t *md, fmd_ituple_t ic_start_send, fmd_ituple_t ic_stop_send,
                       int dim, int dir, ccopier_t ccopy)
{
    fmd_ituple_t ic_src, ic_des;

    LOOP3D(ic_src, ic_start_send, ic_stop_send)
    {
        ic_des[0] = ic_src[0];
        ic_des[1] = ic_src[1];
        ic_des[2] = ic_src[2];

        if (dir < 0)
            ic_des[dim] = ic_src[dim] + md->nc[dim];
        else
            ic_des[dim] = ic_src[dim] - md->nc[dim];

        cell_t *c_src = ARRAY_ELEMENT(md->subd.gridp, ic_src);
        cell_t *c_des = ARRAY_ELEMENT(md->subd.gridp, ic_des);

        if (c_src->parts_num > 0) ccopy(md, c_src, c_des, dim, dir);
    }
}

/* creates an empty buffer */
static void *create_packbuffer_for_partialforces(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop)
{
    fmd_ituple_t ic;
    cell_t *c;
    int s, size = 0;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &s);

    size += s * (ic_stop[0] - ic_start[0]) *
                (ic_stop[1] - ic_start[1]) *
                (ic_stop[2] - ic_start[2]);

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        if (c->parts_num > 0)
        {
            MPI_Pack_size(c->parts_num, md->mpi_types.mpi_rtuple, md->MD_comm, &s);

            size += s;
        }
    }

    return (size > 0 ? m_alloc(md, size) : NULL);
}

static void partialforces_pack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                               int *size, fmd_pointer_t *out, bool nodest)
{
    fmd_ituple_t ic;
    cell_t *c;

    if (nodest) return;

    *out = create_packbuffer_for_partialforces(md, ic_start, ic_stop); /* make an empty buffer */
    *size = 0;

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Pack(&c->parts_num, 1, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);

        if (c->parts_num > 0)
            MPI_Pack(c->F, c->parts_num, md->mpi_types.mpi_rtuple, *out, INT_MAX, size, md->MD_comm);
    }
}

static void partialforces_unpack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                                 int insize, fmd_pointer_t in, int dim, int dir)
{
    fmd_ituple_t ic;
    int byte = 0;

    LOOP3D(ic, ic_start, ic_stop)
    {
        unsigned num;
        cell_t *c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Unpack(in, insize, &byte, &num, 1, MPI_UNSIGNED, md->MD_comm);

        if (num > 0)
        {
            fmd_real_t *F = m_alloc(md, num * sizeof(fmd_rtuple_t));

            MPI_Unpack(in, insize, &byte, F, num, md->mpi_types.mpi_rtuple, md->MD_comm);

            for (int i=0; i < num; i++)
                for (int d=0; d < DIM; d++)
                    FRC(c, i, d) += F[COMP(i, d)];

            free(F);
        }
    }
}

/* creates an empty buffer */
static void *create_packbuffer_for_Femb(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop)
{
    fmd_ituple_t ic;
    cell_t *c;
    int s, size = 0;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &s);

    size += s * (ic_stop[0] - ic_start[0]) *
                (ic_stop[1] - ic_start[1]) *
                (ic_stop[2] - ic_start[2]);

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        if (c->parts_num > 0)
        {
            MPI_Pack_size(c->parts_num, FMD_MPI_REAL, md->MD_comm, &s);

            size += s;
        }
    }

    return (size > 0 ? m_alloc(md, size) : NULL);
}

static void Femb_pack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                      int *size, fmd_pointer_t *out, bool nodest)
{
    fmd_ituple_t ic;
    cell_t *c;

    if (nodest) return;

    *out = create_packbuffer_for_Femb(md, ic_start, ic_stop); /* make an empty buffer */
    *size = 0;

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Pack(&c->parts_num, 1, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);

        if (c->parts_num > 0)
            MPI_Pack(c->vaream, c->parts_num, FMD_MPI_REAL, *out, INT_MAX, size, md->MD_comm);
    }
}

static void Femb_unpack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                        int insize, fmd_pointer_t in, int dim, int dir)
{
    fmd_ituple_t ic;
    int byte = 0;

    LOOP3D(ic, ic_start, ic_stop)
    {
        unsigned num;
        cell_t *c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Unpack(in, insize, &byte, &num, 1, MPI_UNSIGNED, md->MD_comm);

        if (num > 0)
            MPI_Unpack(in, insize, &byte, c->vaream, num, FMD_MPI_REAL, md->MD_comm);
    }
}

/* creates an empty buffer */
static void *create_packbuffer_for_ghostinit(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop)
{
    fmd_ituple_t ic;
    cell_t *c;
    int c_rtuple = 0, c_int = 0, c_unsigned = 0;
    int s, size = 0;

    if (md->cellinfo.x_active) c_rtuple++;
    if (md->cellinfo.GroupID_active) c_int++;
    if (md->cellinfo.atomkind_active) c_unsigned++;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &s);

    size += s * (ic_stop[0] - ic_start[0]) *
                (ic_stop[1] - ic_start[1]) *
                (ic_stop[2] - ic_start[2]);

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        if (c->parts_num > 0)
        {
            if (c_rtuple > 0)
            {
                MPI_Pack_size(c->parts_num, md->mpi_types.mpi_rtuple, md->MD_comm, &s);

                size += c_rtuple * s;
            }

            if (c_int > 0)
            {
                MPI_Pack_size(c->parts_num, MPI_INT, md->MD_comm, &s);

                size += c_int * s;
            }

            if (c_unsigned > 0)
            {
                MPI_Pack_size(c->parts_num, MPI_UNSIGNED, md->MD_comm, &s);

                size += c_unsigned * s;
            }
        }
    }

    return (size > 0 ? m_alloc(md, size) : NULL);
}

static void ghostinit_pack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                           int *size, fmd_pointer_t *out, bool nodest)
{
    fmd_ituple_t ic;
    cell_t *c;

    if (nodest) return;

    *out = create_packbuffer_for_ghostinit(md, ic_start, ic_stop); /* make an empty buffer */
    *size = 0;

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Pack(&c->parts_num, 1, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);

        if (c->parts_num > 0)
        {
            MPI_Pack(c->x, c->parts_num, md->mpi_types.mpi_rtuple, *out, INT_MAX, size, md->MD_comm);
            MPI_Pack(c->GroupID, c->parts_num, MPI_INT, *out, INT_MAX, size, md->MD_comm);
            MPI_Pack(c->atomkind, c->parts_num, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);
        }
    }
}

static inline void init_trans_and_dis(fmd_t *md, int dim, int dir, bool *trans, fmd_real_t *dis)
{
    if (md->subd.is[dim] == 0)
    {
        if (dir == +1)
        {
            *trans = true;
            *dis = -md->l[dim];
        }
    }
    else if (md->subd.is[dim] == md->ns[dim]-1)
    {
        if (dir == -1)
        {
            *trans = true;
            *dis = md->l[dim];
        }
    }
}

static void ghostinit_unpack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                             int insize, fmd_pointer_t in, int dim, int dir)
{
    fmd_ituple_t ic;
    bool trans = false;
    fmd_real_t dis;
    int byte = 0;

    init_trans_and_dis(md, dim, dir, &trans, &dis);

    LOOP3D(ic, ic_start, ic_stop)
    {
        unsigned num;
        cell_t *c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Unpack(in, insize, &byte, &num, 1, MPI_UNSIGNED, md->MD_comm);

        if (num > 0)
        {
            c->parts_num = num;

            if (c->capacity < c->parts_num || c->parts_num + md->cell_inc_dbl < c->capacity)
                _fmd_cell_resize(md, c);

            MPI_Unpack(in, insize, &byte, c->x, num, md->mpi_types.mpi_rtuple, md->MD_comm);

            if (trans)
                for (int i=0; i<num; i++)
                    POS(c, i, dim) += dis;

            MPI_Unpack(in, insize, &byte, c->GroupID, num, MPI_INT, md->MD_comm);

            MPI_Unpack(in, insize, &byte, c->atomkind, num, MPI_UNSIGNED, md->MD_comm);
        }
    }
}

/* creates an empty buffer */
static void *create_packbuffer_for_migrate(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop)
{
    fmd_ituple_t ic;
    cell_t *c;
    int c_rtuple = 0, c_int = 0, c_unsigned = 0;
    int s, size = 0;

    if (md->cellinfo.x_active) c_rtuple++;
    if (md->cellinfo.v_active) c_rtuple++;
    if (md->cellinfo.GroupID_active) c_int++;
    if (md->cellinfo.AtomID_active) c_unsigned++;
    if (md->cellinfo.atomkind_active) c_unsigned++;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &s);

    size += s * (ic_stop[0] - ic_start[0]) *
                (ic_stop[1] - ic_start[1]) *
                (ic_stop[2] - ic_start[2]);

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        if (c->parts_num > 0)
        {
            if (c_rtuple > 0)
            {
                MPI_Pack_size(c->parts_num, md->mpi_types.mpi_rtuple, md->MD_comm, &s);

                size += c_rtuple * s;
            }

            if (c_int > 0)
            {
                MPI_Pack_size(c->parts_num, MPI_INT, md->MD_comm, &s);

                size += c_int * s;
            }

            if (c_unsigned > 0)
            {
                MPI_Pack_size(c->parts_num, MPI_UNSIGNED, md->MD_comm, &s);

                size += c_unsigned * s;
            }
        }
    }

    return (size > 0 ? m_alloc(md, size) : NULL);
}

static void migrate_pack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                         int *size, fmd_pointer_t *out, bool nodest)
{
    fmd_ituple_t ic;
    cell_t *c;

    if (nodest) /* no destination; simply remove the atoms */
    {
        md->subd.NumberOfParticles -= set_partsnum_to_zero(md, ic_start, ic_stop);

        return;
    }

    *out = create_packbuffer_for_migrate(md, ic_start, ic_stop); /* make an empty buffer */
    *size = 0;

    LOOP3D(ic, ic_start, ic_stop)
    {
        c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Pack(&c->parts_num, 1, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);

        if (c->parts_num > 0)
        {
            MPI_Pack(c->x, c->parts_num, md->mpi_types.mpi_rtuple, *out, INT_MAX, size, md->MD_comm);

            if (c->v != NULL)
                MPI_Pack(c->v, c->parts_num, md->mpi_types.mpi_rtuple, *out, INT_MAX, size, md->MD_comm);

            if (c->GroupID != NULL)
                MPI_Pack(c->GroupID, c->parts_num, MPI_INT, *out, INT_MAX, size, md->MD_comm);

            if (c->AtomID != NULL)
                MPI_Pack(c->AtomID, c->parts_num, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);

            MPI_Pack(c->atomkind, c->parts_num, MPI_UNSIGNED, *out, INT_MAX, size, md->MD_comm);

            md->subd.NumberOfParticles -= c->parts_num;
            c->parts_num = 0;
        }
    }
}

static void migrate_unpack(fmd_t *md, fmd_ituple_t ic_start, fmd_ituple_t ic_stop,
                           int insize, fmd_pointer_t in, int dim, int dir)
{
    fmd_ituple_t ic;
    bool trans = false;
    fmd_real_t dis;
    int byte = 0;

    init_trans_and_dis(md, dim, dir, &trans, &dis);

    LOOP3D(ic, ic_start, ic_stop)
    {
        unsigned incr;
        cell_t *c = ARRAY_ELEMENT(md->subd.gridp, ic);

        MPI_Unpack(in, insize, &byte, &incr, 1, MPI_UNSIGNED, md->MD_comm);

        if (incr > 0)
        {
            unsigned oldnum = c->parts_num;

            c->parts_num += incr;
            if (c->parts_num > c->capacity) _fmd_cell_resize(md, c);

            md->subd.NumberOfParticles += incr;

            MPI_Unpack(in, insize, &byte, c->x + oldnum*DIM, incr, md->mpi_types.mpi_rtuple, md->MD_comm);

            if (trans)
                for (int i=0; i<incr; i++)
                    POS(c, oldnum + i, dim) += dis;

            if (c->v != NULL)
                MPI_Unpack(in, insize, &byte, c->v + oldnum*DIM, incr, md->mpi_types.mpi_rtuple, md->MD_comm);

            if (c->GroupID != NULL)
                MPI_Unpack(in, insize, &byte, c->GroupID + oldnum, incr, MPI_INT, md->MD_comm);

            if (c->AtomID != NULL)
                MPI_Unpack(in, insize, &byte, c->AtomID + oldnum, incr, MPI_UNSIGNED, md->MD_comm);

            MPI_Unpack(in, insize, &byte, c->atomkind + oldnum, incr, MPI_UNSIGNED, md->MD_comm);
        }
    }
}

static void transfer(fmd_t *md, int source, int dest,
                     fmd_ituple_t ic_start_send, fmd_ituple_t ic_stop_send,
                     fmd_ituple_t ic_start_recv, fmd_ituple_t ic_stop_recv,
                     packer_t pack, unpacker_t unpack, int dim, int dir)
{
    int send_size, recv_size;
    void *send_buff, *recv_buff;
    MPI_Request request;
    MPI_Status status;

    pack(md, ic_start_send, ic_stop_send, &send_size, &send_buff, dest == MPI_PROC_NULL);

    if (dest != MPI_PROC_NULL) MPI_Isend(send_buff, send_size, MPI_PACKED, dest, 2, md->MD_comm, &request);

    if (source != MPI_PROC_NULL)
    {
        MPI_Probe(source, 2, md->MD_comm, &status);
        MPI_Get_count(&status, MPI_PACKED, &recv_size);
        recv_buff = m_alloc(md, recv_size);
        MPI_Recv(recv_buff, recv_size, MPI_PACKED, source, 2, md->MD_comm, &status);

        unpack(md, ic_start_recv, ic_stop_recv, recv_size, recv_buff, dim, dir);

        free(recv_buff);
    }

    if (dest != MPI_PROC_NULL)
    {
        MPI_Wait(&request, &status);

        free(send_buff);
    }
}

static void set_cell_indexes_for_migration_in_direction_d(
    Subdomain_t *s, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    int dd;

    for (dd=0; dd<DIM; dd++)
    {
        if (dd == d) // only ghost
        {
            ic_start_send_lower[dd] = s->ic_start[dd];
            ic_stop_send_lower[dd] = ic_start_send_lower[dd] + s->ic_start[dd];
            ic_start_receive_lower[dd] = 0;
            ic_stop_receive_lower[dd] = s->ic_start[dd];
            ic_stop_send_upper[dd] = s->ic_stop[dd];
            ic_start_send_upper[dd] = ic_stop_send_upper[dd] - s->ic_start[dd];
            ic_start_receive_upper[dd] = s->ic_stop[dd];
            ic_stop_receive_upper[dd] = s->cell_num[dd];
        }
        else if (dd > d) // including ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd] = 0;
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd] = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd] = s->cell_num[dd];
        }
        else // excluding ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd] = s->ic_start[dd];
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd] = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd] = s->ic_stop[dd];
        }
    }
}

static void set_cell_indexes_for_transfer_in_direction_d(
    Subdomain_t *s, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    int dd;

    for (dd=0; dd<DIM; dd++)
    {
        if (dd == d) // only ghost
        {
            ic_start_send_lower[dd] = s->ic_start[dd];
            ic_stop_send_lower[dd] = ic_start_send_lower[dd] + 1;
            ic_stop_receive_lower[dd] = s->ic_start[dd];
            ic_start_receive_lower[dd] = ic_stop_receive_lower[dd] - 1;
            ic_stop_send_upper[dd] = s->ic_stop[dd];
            ic_start_send_upper[dd] = ic_stop_send_upper[dd] - 1;
            ic_start_receive_upper[dd] = s->ic_stop[dd];
            ic_stop_receive_upper[dd] = ic_start_receive_upper[dd] + 1;
        }
        else if (dd > d) // including ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd]
                                       = s->ic_start[dd] - 1;
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd]
                                       = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd]
                                       = s->ic_stop[dd] + 1;
        }
        else // excluding ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd]
                                       = s->ic_start[dd];
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd]
                                       = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd]
                                       = s->ic_stop[dd];
        }
    }
}

void _fmd_ghostparticles_transfer_Femb(fmd_t *md)
{
    int d;
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (d = DIM-1; d >= 0; d--)
    {
        set_cell_indexes_for_transfer_in_direction_d(
            &md->subd, d, ic_start_send_lower, ic_stop_send_lower,
            ic_start_receive_lower, ic_stop_receive_lower,
            ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
            ic_stop_receive_upper);

        if (md->ns[d] != 1)
        {
            /* sending to lower process, receiving from upper process */

            transfer(md, md->subd.rank_of_upper_subd[d], md->subd.rank_of_lower_subd[d],
                     ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper, ic_stop_receive_upper,
                     Femb_pack, Femb_unpack, d, -1);

            /* sending to upper process, receiving from lower process */

            transfer(md, md->subd.rank_of_lower_subd[d], md->subd.rank_of_upper_subd[d],
                     ic_start_send_upper, ic_stop_send_upper, ic_start_receive_lower, ic_stop_receive_lower,
                     Femb_pack, Femb_unpack, d, +1);
        }
        else if (md->PBC[d])
        {
            copy_local(md, ic_start_send_lower, ic_stop_send_lower, d, -1, ccopy_for_Femb);
            copy_local(md, ic_start_send_upper, ic_stop_send_upper, d, +1, ccopy_for_Femb);
        }

    }
}

void _fmd_ghostparticles_init(fmd_t *md)
{
    int d;
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (d = DIM-1; d >= 0; d--)
    {
        set_cell_indexes_for_transfer_in_direction_d(
            &md->subd, d, ic_start_send_lower, ic_stop_send_lower,
            ic_start_receive_lower, ic_stop_receive_lower,
            ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
            ic_stop_receive_upper);

        if (md->ns[d] != 1)
        {
            /* sending to lower process, receiving from upper process */

            transfer(md, md->subd.rank_of_upper_subd[d], md->subd.rank_of_lower_subd[d],
                     ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper, ic_stop_receive_upper,
                     ghostinit_pack, ghostinit_unpack, d, -1);

            /* sending to upper process, receiving from lower process */

            transfer(md, md->subd.rank_of_lower_subd[d], md->subd.rank_of_upper_subd[d],
                     ic_start_send_upper, ic_stop_send_upper, ic_start_receive_lower, ic_stop_receive_lower,
                     ghostinit_pack, ghostinit_unpack, d, +1);
        }
        else if (md->PBC[d])
        {
            copy_local(md, ic_start_send_lower, ic_stop_send_lower, d, -1, ccopy_for_ghostinit);
            copy_local(md, ic_start_send_upper, ic_stop_send_upper, d, +1, ccopy_for_ghostinit);
        }
    }
}

void _fmd_ghostparticles_transfer_partialforces(fmd_t *md)
{
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (int d = 0; d < DIM; d++)
    {
        set_cell_indexes_for_transfer_in_direction_d(
            &md->subd, d, ic_start_receive_lower, ic_stop_receive_lower,
            ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper,
            ic_stop_receive_upper, ic_start_send_upper, ic_stop_send_upper );

        if (md->ns[d] != 1)
        {
            /* sending to lower process, receiving from upper process */

            transfer(md, md->subd.rank_of_upper_subd[d], md->subd.rank_of_lower_subd[d],
                     ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper, ic_stop_receive_upper,
                     partialforces_pack, partialforces_unpack, d, -1);

            /* sending to upper process, receiving from lower process */

            transfer(md, md->subd.rank_of_lower_subd[d], md->subd.rank_of_upper_subd[d],
                     ic_start_send_upper, ic_stop_send_upper, ic_start_receive_lower, ic_stop_receive_lower,
                     partialforces_pack, partialforces_unpack, d, +1);
        }
        else
        {
            if (md->PBC[d])
            {
                copy_local(md, ic_start_send_lower, ic_stop_send_lower, d, -1, ccopy_for_partialforces);
                copy_local(md, ic_start_send_upper, ic_stop_send_upper, d, +1, ccopy_for_partialforces);
            }
        }
    }
}

void _fmd_particles_migrate(fmd_t *md)
{
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (int d = 0; d < DIM; d++)
    {
        set_cell_indexes_for_migration_in_direction_d(
            &md->subd, d, ic_start_receive_lower, ic_stop_receive_lower,
            ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper,
            ic_stop_receive_upper, ic_start_send_upper, ic_stop_send_upper );

        if (md->ns[d] != 1)
        {
            /* sending to lower process, receiving from upper process */

            transfer(md, md->subd.rank_of_upper_subd[d], md->subd.rank_of_lower_subd[d],
                     ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper, ic_stop_receive_upper,
                     migrate_pack, migrate_unpack, d, -1);

            /* sending to upper process, receiving from lower process */

            transfer(md, md->subd.rank_of_lower_subd[d], md->subd.rank_of_upper_subd[d],
                     ic_start_send_upper, ic_stop_send_upper, ic_start_receive_lower, ic_stop_receive_lower,
                     migrate_pack, migrate_unpack, d, +1);
        }
        else
        {
            if (md->PBC[d])
            {
                copy_local(md, ic_start_send_lower, ic_stop_send_lower, d, -1, ccopy_for_migrate);
                copy_local(md, ic_start_send_upper, ic_stop_send_upper, d, +1, ccopy_for_migrate);
            }
            else
            {
                md->subd.NumberOfParticles -= set_partsnum_to_zero(md, ic_start_send_lower, ic_stop_send_lower);
                md->subd.NumberOfParticles -= set_partsnum_to_zero(md, ic_start_send_upper, ic_stop_send_upper);
            }
        }
    }
}
