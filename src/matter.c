/*
  matter.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2022 Arham Amouye Foumani

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

#include <tgmath.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include "matter.h"
#include "fmd-private.h"
#include "general.h"
#include "molecule.h"

/* calculate GroupTemperature, GroupParticlesNum and GroupKineticEnergy from "local grid" */
void _fmd_compute_GroupTemperature_etc_localgrid(fmd_t *md)
{
    fmd_ituple_t ic;
    fmd_real_t m_vSqd_Sum = 0, m_vSqd_SumSum;
    fmd_real_t mass;
    int ParticlesNum = 0;
    fmd_rtuple_t MomentumSum = {0., 0., 0.};
    cell_t *c;
    unsigned i;

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            ParticlesNum++;

            mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            for (int d=0; d<DIM; d++)
                MomentumSum[d] += mass * VEL(c, i, d);

            m_vSqd_Sum += mass * ( sqrr(VEL(c, i, 0)) +
                                   sqrr(VEL(c, i, 1)) +
                                   sqrr(VEL(c, i, 2)) );
        }

    MPI_Allreduce(MomentumSum, md->GroupMomentum, DIM, FMD_MPI_REAL, MPI_SUM, md->MD_comm);

    MPI_Allreduce(&ParticlesNum, &md->GroupParticlesNum, 1, MPI_INT, MPI_SUM, md->MD_comm);

    MPI_Allreduce(&m_vSqd_Sum, &m_vSqd_SumSum, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
    md->GroupKineticEnergy = 0.5 * m_vSqd_SumSum;

    md->GroupTemperature = m_vSqd_SumSum / (3.0 * md->GroupParticlesNum * K_BOLTZMANN);
}

/* calculate GroupTemperature, GroupParticlesNum and GroupKineticEnergy from "global grid" */
static void calculate_GroupTemperature_etc(fmd_t *md)
{
    if (md->Is_MD_comm_root)
    {
        fmd_ituple_t ic;
        fmd_real_t m_vSqd_Sum = 0.0;
        cell_t *c;
        int i;

        md->GroupParticlesNum = 0;
        md->GroupMomentum[0] = md->GroupMomentum[1] = md->GroupMomentum[2] = 0.0;

        LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
            for (c=&ARRAY_ELEMENT(md->global_grid, ic), i=0; i < c->parts_num; i++)
            {
                if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                    continue;

                md->GroupParticlesNum++;

                fmd_real_t mass = md->potsys.atomkinds[c->atomkind[i]].mass;

                for (int d=0; d<DIM; d++)
                {
                    m_vSqd_Sum += mass * sqrr(VEL(c, i, d));
                    md->GroupMomentum[d] += mass * VEL(c, i, d);
                }
            }

        md->GroupTemperature = m_vSqd_Sum / (3.0 * md->GroupParticlesNum * K_BOLTZMANN);
    }

    MPI_Bcast(&md->GroupParticlesNum, 1, MPI_UNSIGNED, RANK0, md->MD_comm);
    MPI_Bcast(&md->GroupTemperature, 1, FMD_MPI_REAL, RANK0, md->MD_comm);
    MPI_Bcast(md->GroupMomentum, DIM, FMD_MPI_REAL, RANK0, md->MD_comm);

    md->GroupKineticEnergy = 3.0/2.0 * md->GroupParticlesNum * K_BOLTZMANN * md->GroupTemperature;
}

void fmd_matt_addVelocity(fmd_t *md, int GroupID, fmd_real_t vx, fmd_real_t vy, fmd_real_t vz)
{
    cell_t ***grid;
    const int *start, *stop;
    fmd_ituple_t ic;

    if (md->ParticlesDistributed)
    {
        grid = md->Subdomain.grid;
        start = md->Subdomain.ic_start;
        stop = md->Subdomain.ic_stop;
    }
    else
    {
        start = _fmd_ThreeZeros_int;
        if (md->Is_MD_comm_root)
        {
            grid = md->global_grid;
            stop = md->nc;
        }
        else
            stop = start;
    }

    int i;
    cell_t *c;

    LOOP3D(ic, start, stop)
        for (c = &ARRAY_ELEMENT(grid, ic), i=0; i < c->parts_num; i++)
            if (GroupID == FMD_GROUP_ALL || GroupID == c->GroupID[i])
            {
                VEL(c, i, 0) += vx;
                VEL(c, i, 1) += vy;
                VEL(c, i, 2) += vz;
            }

    if (md->ParticlesDistributed)
        _fmd_compute_GroupTemperature_etc_localgrid(md);
    else
        calculate_GroupTemperature_etc(md);
}

void fmd_matt_findLimits(fmd_t *md, fmd_rtuple_t LowerLimit, fmd_rtuple_t UpperLimit)
{
    cell_t ***grid;
    const int *start, *stop;

    if (md->ParticlesDistributed)
    {
        grid = md->Subdomain.grid;
        start = md->Subdomain.ic_start;
        stop = md->Subdomain.ic_stop;
    }
    else
    {
        start = _fmd_ThreeZeros_int;
        if (md->Is_MD_comm_root)
        {
            grid = md->global_grid;
            stop = md->nc;
        }
        else
            stop = start;
    }

    fmd_rtuple_t L, U;

    for (int d=0; d<DIM; d++)
    {
        L[d] = DBL_MAX;
        U[d] = DBL_MIN;
    }

    int i;
    cell_t *c;
    fmd_ituple_t ic;

    LOOP3D(ic, start, stop)
        for (c = &ARRAY_ELEMENT(grid, ic), i=0; i < c->parts_num; i++)
            for (int d=0; d<DIM; d++)
            {
                if (POS(c, i, d) < L[d])
                    L[d] = POS(c, i, d);
                if (POS(c, i, d) > U[d])
                    U[d] = POS(c, i, d);
            }

    if (md->ParticlesDistributed)
    {
        MPI_Allreduce(L, LowerLimit, DIM, FMD_MPI_REAL, MPI_MIN, md->MD_comm);
        MPI_Allreduce(U, UpperLimit, DIM, FMD_MPI_REAL, MPI_MAX, md->MD_comm);
    }
    else
    {
        MPI_Bcast(L, DIM, FMD_MPI_REAL, RANK0, md->MD_comm);
        MPI_Bcast(U, DIM, FMD_MPI_REAL, RANK0, md->MD_comm);

        for (int d=0; d<DIM; d++)
        {
            LowerLimit[d] = L[d];
            UpperLimit[d] = U[d];
        }
    }
}

static void find_global_start_stop_ic_of_a_subd(fmd_t *md, int rank,
 fmd_ituple_t global_icstart, fmd_ituple_t global_icstop)
{
    fmd_ituple_t is;
    int r, w;

    INDEX_3D(rank, md->ns, is);

    for (int d=0; d<3; d++)
    {
        r = md->nc[d] % md->ns[d];
        w = md->nc[d] / md->ns[d];

        if (is[d] < r)
        {
            global_icstart[d] = is[d] * (w + 1);
            global_icstop[d] = global_icstart[d] + w + 1;
        }
        else
        {
            global_icstart[d] = is[d] * w + r;
            global_icstop[d] = global_icstart[d] + w;
        }
    }
}

/* creates an empty buffer */
static void *create_packbuffer_for_matt_distribute(fmd_t *md, fmd_ituple_t global_icstart, fmd_ituple_t global_icstop)
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
    if (md->cellinfo.molkind_active) c_unsigned += 3; /* for molkind, MolID and AtomIDlocal arrays */

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &s);

    size += s * (global_icstop[0] - global_icstart[0]) *
                (global_icstop[1] - global_icstart[1]) *
                (global_icstop[2] - global_icstart[2]);

    LOOP3D(ic, global_icstart, global_icstop)
    {
        c = &ARRAY_ELEMENT(md->global_grid, ic);

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

    return (size > 0 ? m_alloc(size) : NULL);
}

static void pack_for_matt_distribute(fmd_t *md, void *buff, int *bytecount,
                                     fmd_ituple_t global_icstart, fmd_ituple_t global_icstop)
{
    *bytecount = 0;

    fmd_ituple_t ic;
    cell_t *c;

    LOOP3D(ic, global_icstart, global_icstop)
    {
        c = &ARRAY_ELEMENT(md->global_grid, ic);

        MPI_Pack(&c->parts_num, 1, MPI_UNSIGNED, buff, INT_MAX, bytecount, md->MD_comm);

        if (c->parts_num > 0)
        {
            if (c->x != NULL)
                MPI_Pack(c->x, c->parts_num, md->mpi_types.mpi_rtuple, buff, INT_MAX, bytecount, md->MD_comm);

            if (c->v != NULL)
                MPI_Pack(c->v, c->parts_num, md->mpi_types.mpi_rtuple, buff, INT_MAX, bytecount, md->MD_comm);

            if (c->GroupID != NULL)
                MPI_Pack(c->GroupID, c->parts_num, MPI_INT, buff, INT_MAX, bytecount, md->MD_comm);

            if (c->AtomID != NULL)
                MPI_Pack(c->AtomID, c->parts_num, MPI_UNSIGNED, buff, INT_MAX, bytecount, md->MD_comm);

            if (c->atomkind != NULL)
                MPI_Pack(c->atomkind, c->parts_num, MPI_UNSIGNED, buff, INT_MAX, bytecount, md->MD_comm);

            if (c->molkind != NULL)
            {
                MPI_Pack(c->molkind, c->parts_num, MPI_UNSIGNED, buff, INT_MAX, bytecount, md->MD_comm);
                MPI_Pack(c->MolID, c->parts_num, MPI_UNSIGNED, buff, INT_MAX, bytecount, md->MD_comm);
                MPI_Pack(c->AtomIDlocal, c->parts_num, MPI_UNSIGNED, buff, INT_MAX, bytecount, md->MD_comm);
            }
        }

        _fmd_cell_free(c);
    }
}

static void transfer_from_globalgrid_to_rank0_grid(fmd_t *md)
{
    fmd_ituple_t icl, icg;
    cell_t *cl, *cg;

    LOOP3D(icl, md->Subdomain.ic_start, md->Subdomain.ic_stop)
    {
        _fmd_conv_ic_loc_to_glob(md, icl, icg);

        cl = &ARRAY_ELEMENT(md->Subdomain.grid, icl);
        cg = &ARRAY_ELEMENT(md->global_grid, icg);

        cl->parts_num = cg->parts_num;
        md->Subdomain.NumberOfParticles += cl->parts_num;

        _fmd_cell_resize(md, cl);

        if (cg->x != NULL)
            memcpy(cl->x, cg->x, cg->parts_num * sizeof(fmd_rtuple_t));

        if (cg->v != NULL)
            memcpy(cl->v, cg->v, cg->parts_num * sizeof(fmd_rtuple_t));

        if (cg->GroupID != NULL)
            memcpy(cl->GroupID, cg->GroupID, cg->parts_num * sizeof(int));

        if (cg->AtomID != NULL)
            memcpy(cl->AtomID, cg->AtomID, cg->parts_num * sizeof(unsigned));

        if (cg->atomkind != NULL)
            memcpy(cl->atomkind, cg->atomkind, cg->parts_num * sizeof(unsigned));

        if (cg->molkind != NULL)
        {
            memcpy(cl->molkind, cg->molkind, cg->parts_num * sizeof(unsigned));
            memcpy(cl->MolID, cg->MolID, cg->parts_num * sizeof(unsigned));
            memcpy(cl->AtomIDlocal, cg->AtomIDlocal, cg->parts_num * sizeof(unsigned));
        }

        _fmd_cell_free(cg);
    }
}

static void unpack_for_matt_distribute(fmd_t *md, void *packbuf, int bufsize)
{
    fmd_ituple_t ic;
    cell_t *c;
    int pos = 0;

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
    {
        c = &ARRAY_ELEMENT(md->Subdomain.grid, ic);

        MPI_Unpack(packbuf, bufsize, &pos, &c->parts_num, 1, MPI_UNSIGNED, md->MD_comm);

        if (c->parts_num > 0)
        {
            md->Subdomain.NumberOfParticles += c->parts_num;

            _fmd_cell_resize(md, c);

            if (c->x != NULL)
                MPI_Unpack(packbuf, bufsize, &pos, c->x, c->parts_num, md->mpi_types.mpi_rtuple, md->MD_comm);

            if (c->v != NULL)
                MPI_Unpack(packbuf, bufsize, &pos, c->v, c->parts_num, md->mpi_types.mpi_rtuple, md->MD_comm);

            if (c->GroupID != NULL)
                MPI_Unpack(packbuf, bufsize, &pos, c->GroupID, c->parts_num, MPI_INT, md->MD_comm);

            if (c->AtomID != NULL)
                MPI_Unpack(packbuf, bufsize, &pos, c->AtomID, c->parts_num, MPI_UNSIGNED, md->MD_comm);

            if (c->atomkind != NULL)
                MPI_Unpack(packbuf, bufsize, &pos, c->atomkind, c->parts_num, MPI_UNSIGNED, md->MD_comm);

            if (c->molkind != NULL)
            {
                MPI_Unpack(packbuf, bufsize, &pos, c->molkind, c->parts_num, MPI_UNSIGNED, md->MD_comm);
                MPI_Unpack(packbuf, bufsize, &pos, c->MolID, c->parts_num, MPI_UNSIGNED, md->MD_comm);
                MPI_Unpack(packbuf, bufsize, &pos, c->AtomIDlocal, c->parts_num, MPI_UNSIGNED, md->MD_comm);
            }
        }
    }
}

void _fmd_matt_distribute(fmd_t *md)
{
    if (md->Subdomain.grid == NULL) _fmd_subd_init(md);

    calculate_GroupTemperature_etc(md);

    if (md->Is_MD_comm_root)
    {
        for (int i=1; i < md->Subdomain.numprocs; i++)  /* for all processes execpt the root */
        {
            fmd_ituple_t global_icstart, global_icstop;
            int bytecount;

            find_global_start_stop_ic_of_a_subd(md, i, global_icstart, global_icstop);

            void *buff = create_packbuffer_for_matt_distribute(md, global_icstart, global_icstop);

            pack_for_matt_distribute(md, buff, &bytecount, global_icstart, global_icstop);

            MPI_Send(buff, bytecount, MPI_PACKED, i, 51, md->MD_comm);

            free(buff);
        }

        transfer_from_globalgrid_to_rank0_grid(md); /* for the root process */

        _fmd_array_ordinary3d_free((void ***)md->global_grid, md->nc[0], md->nc[1]);
    }
    else
    {
        MPI_Status status;
        int bufsize;

        MPI_Probe(RANK0, 51, md->MD_comm, &status);

        MPI_Get_count(&status, MPI_PACKED, &bufsize);

        void *buff = m_alloc(bufsize);

        MPI_Recv(buff, bufsize, MPI_PACKED, RANK0, 51, md->MD_comm, &status);

        unpack_for_matt_distribute(md, buff, bufsize);

        free(buff);
    }

    MPI_Bcast(&md->TotalNoOfParticles, 1, MPI_UNSIGNED, RANK0, md->MD_comm);
    MPI_Bcast(&md->TotalNoOfMolecules, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

    if (md->TotalNoOfMolecules > 0) _fmd_matt_updateAtomNeighbors(md);

    md->GlobalGridExists = false;
    md->ParticlesDistributed = true;
}

void fmd_matt_changeGroupID(fmd_t *md, int old, int new)
{
    int i;
    cell_t *c;
    fmd_ituple_t ic;

    if (md->ParticlesDistributed)
    {
        LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
            for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < c->parts_num; i++)
                if (old == FMD_GROUP_ALL || c->GroupID[i] == old) c->GroupID[i] = new;
    }
    else
    {
        if (md->Is_MD_comm_root)
        {
            LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
                for (c = &ARRAY_ELEMENT(md->global_grid, ic), i=0; i < c->parts_num; i++)
                    if (old == FMD_GROUP_ALL || c->GroupID[i] == old) c->GroupID[i] = new;
        }
    }
}

#define WHAT_IF_THE_PARTICLE_HAS_LEFT_THE_BOX_2(md, c, i, d, x)                \
    if ( (((x) < 0.0) || ((x) >= (md)->l[(d)])) )                              \
    {                                                                          \
        if (!(md)->PBC[(d)])                                                   \
        {                                                                      \
            _fmd_cell_remove_atom((md), (c), (i));                             \
            (md)->Subdomain.NumberOfParticles--;                               \
            (i)--;                                                             \
            break;                                                             \
        }                                                                      \
        else                                                                   \
            (x) += ((x) < 0.0 ? (md)->l[(d)] : -(md)->l[(d)]);                 \
    }                                                                          \
    do {} while (0)

void fmd_matt_translate(fmd_t *md, int GroupID, fmd_real_t dx, fmd_real_t dy, fmd_real_t dz)
{
    assert(!md->ParticlesDistributed);

    if (!md->Is_MD_comm_root) return;

    fmd_rtriple_t dr;

    dr[0] = dx;
    dr[1] = dy;
    dr[2] = dz;

    fmd_ituple_t ic, jc;
    int i;
    cell_t *c;

    LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
        for (c = &ARRAY_ELEMENT(md->global_grid, ic), i=0; i < c->parts_num; i++)
            if (GroupID == FMD_GROUP_ALL || GroupID == c->GroupID[i])
                for (int d=0; d<3; d++)
                {
                    POS(c, i, d) += dr[d];

                    WHAT_IF_THE_PARTICLE_HAS_LEFT_THE_BOX_2(md, c, i, d, POS(c, i, d));
                }

    LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
        for (c = &ARRAY_ELEMENT(md->global_grid, ic), i=0; i < c->parts_num; i++)
            if (GroupID == FMD_GROUP_ALL || GroupID == c->GroupID[i])
            {
                for (int d=0; d<3; d++)
                {
                    jc[d] = (int)floor(POS(c, i, d) / md->cellh[d]);

                    assert(!(jc[d] < 0 || jc[d] >= md->nc[d])); /* TO-DO: handle error */
                }

                if ((ic[0] != jc[0]) || (ic[1] != jc[1]) || (ic[2] != jc[2]))
                {
                    cell_t *c2 = &ARRAY_ELEMENT(md->global_grid, jc);
                    unsigned j = _fmd_cell_new_particle(md, c2);

                    _fmd_cell_copy_atom_from_cell_to_cell(c, i, c2, j);
                    _fmd_cell_remove_atom(md, c, i);

                    i--;
                }
            }
}

fmd_real_t fmd_matt_getTotalEnergy(fmd_t *md)
{
    return md->GroupKineticEnergy + md->GroupPotentialEnergy;
}

void fmd_matt_giveTemperature(fmd_t *md, int GroupID)
{
    cell_t ***grid;
    const int *start, *stop;
    fmd_ituple_t ic;

    if (md->ParticlesDistributed)
    {
        grid = md->Subdomain.grid;
        start = md->Subdomain.ic_start;
        stop = md->Subdomain.ic_stop;
    }
    else
    {
        start = _fmd_ThreeZeros_int;
        if (md->Is_MD_comm_root)
        {
            grid = md->global_grid;
            stop = md->nc;
        }
        else
            stop = start;
    }

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    cell_t *c;
    int pi;

    LOOP3D(ic, start, stop)
        for (c = &ARRAY_ELEMENT(grid, ic), pi = 0; pi < c->parts_num; pi++)
            if (GroupID == FMD_GROUP_ALL || GroupID == c->GroupID[pi])
            {
                fmd_real_t mass = md->potsys.atomkinds[c->atomkind[pi]].mass;
                fmd_real_t StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

                for (int d=0; d<DIM; d++)
                    VEL(c, pi, d) = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
            }

    gsl_rng_free(rng);
}

fmd_real_t fmd_matt_getTemperature(fmd_t *md)
{
    return md->GroupTemperature;
}

void fmd_matt_getMomentum(fmd_t *md, fmd_rtuple_t out)
{
    for (int d=0; d<DIM; d++)
        out[d] = md->GroupMomentum[d];
}

static config_atom_t *prepare_localdata_for_saveconfig(fmd_t *md)
{
    config_atom_t *localdata;
    fmd_ituple_t ic;
    cell_t *c;
    int pind;
    int k = 0;

    localdata = (config_atom_t *)m_alloc(md->Subdomain.NumberOfParticles * sizeof(config_atom_t));

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), pind=0; pind < c->parts_num; pind++)
        {
            localdata[k].x[0] = (float)POS(c, pind, 0);
            localdata[k].x[1] = (float)POS(c, pind, 1);
            localdata[k].x[2] = (float)POS(c, pind, 2);

            localdata[k].var = c->GroupID[pind];
            localdata[k].atomkind = c->atomkind[pind];

            k++;
        }

    return localdata;
}

static config_atom_t *gather_localdata_on_root_for_saveconfig(fmd_t *md, config_atom_t *localdata)
{
    int *displs;
    unsigned *nums;
    config_atom_t *globaldata;

    if (md->Is_MD_comm_root) nums = (unsigned *)m_alloc(md->Subdomain.numprocs * sizeof(unsigned));
    MPI_Gather(&md->Subdomain.NumberOfParticles, 1, MPI_UNSIGNED, nums, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        md->TotalNoOfParticles = 0;

        for (int i=0; i < md->Subdomain.numprocs; i++)
            md->TotalNoOfParticles += nums[i];

        int displ = 0;

        globaldata = (config_atom_t *)m_alloc(md->TotalNoOfParticles * sizeof(config_atom_t));
        displs     = (int *)m_alloc(md->Subdomain.numprocs * sizeof(int));

        for (int k=0; k < md->Subdomain.numprocs; k++)
        {
            displs[k] = displ;
            displ += nums[k];
        }
    }

    MPI_Gatherv(localdata, md->Subdomain.NumberOfParticles, md->mpi_types.mpi_configa,
        globaldata, nums, displs, md->mpi_types.mpi_configa, RANK0, md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        free(nums);
        free(displs);
    }

    return globaldata;
}

static void save_XYZ_data(fmd_t *md, config_atom_t *globaldata)
{
    char *ElementName;

    fprintf(md->ConfigFilep, "%d\n\n", md->TotalNoOfParticles);

    for (int i=0; i < md->potsys.atomkinds_num; i++)
    {
        ElementName = md->potsys.atomkinds[i].name;

        for (int k=0; k < md->TotalNoOfParticles; k++)
            if (globaldata[k].atomkind == i)
                fprintf(md->ConfigFilep, "%s\t%.2f\t%.2f\t%.2f\n", ElementName,
                        globaldata[k].x[0], globaldata[k].x[1], globaldata[k].x[2]);
    }
}

static void save_VTF_file(fmd_t *md, config_atom_t *globaldata)
{
    char ConfigPath[MAX_PATH_LENGTH];
    int AtomID = 0;

    sprintf(ConfigPath, "%s%05d.vtf", md->SaveDirectory, md->_FileIndex++);
    md->ConfigFilep = f_open(ConfigPath, "w");

    for (int i=0; i < md->potsys.atomkinds_num; i++)
    {
        int count = 0;
        for (int k=0; k < md->TotalNoOfParticles; k++)
            if (globaldata[k].atomkind == i)
                count++;
        fprintf(md->ConfigFilep, "atom %d:%d\tname %s\n", AtomID,
                AtomID+count-1, md->potsys.atomkinds[i].name);
        AtomID += count;
    }

    AtomID = 0;
    for (int i=0; i < md->potsys.atomkinds_num; i++)
    {
        for (int k=0; k < md->TotalNoOfParticles; k++)
            if (globaldata[k].atomkind == i)
            {
                fprintf(md->ConfigFilep, "%d\tbeta %.4f\n", AtomID, globaldata[k].var);
                AtomID++;
            }
    }

    fprintf(md->ConfigFilep, "timestep\n");

    for (int i=0; i < md->potsys.atomkinds_num; i++)
    {
        for (int k=0; k < md->TotalNoOfParticles; k++)
            if (globaldata[k].atomkind == i)
                fprintf(md->ConfigFilep, "%.2f\t%.2f\t%.2f\n",
                    globaldata[k].x[0], globaldata[k].x[1], globaldata[k].x[2]);
    }

    fclose(md->ConfigFilep);
}

void fmd_matt_saveConfiguration(fmd_t *md)
{
    config_atom_t *localdata, *globaldata;

    localdata = prepare_localdata_for_saveconfig(md);
    globaldata = gather_localdata_on_root_for_saveconfig(md, localdata);

    free(localdata);

    if (!md->Is_MD_comm_root) return;

    char ConfigPath[MAX_PATH_LENGTH];

    switch (md->SaveConfigMode)
    {
        case FMD_SCM_XYZ_ATOMSNUM:
            if (md->TotalNoOfParticles != md->_OldNumberOfParticles)
            {
                if (md->_OldNumberOfParticles != -1) fclose(md->ConfigFilep);
                sprintf(ConfigPath, "%s%d.xyz", md->SaveDirectory, md->TotalNoOfParticles);
                md->ConfigFilep = f_open(ConfigPath, "w");
                md->_OldNumberOfParticles = md->TotalNoOfParticles;
            }
            save_XYZ_data(md, globaldata);
            fflush(md->ConfigFilep);
            break;

        case FMD_SCM_XYZ_SEPARATE:
            sprintf(ConfigPath, "%s%05d.xyz", md->SaveDirectory, md->_FileIndex++);
            md->ConfigFilep = f_open(ConfigPath, "w");
            save_XYZ_data(md, globaldata);
            fclose(md->ConfigFilep);
            break;

        case FMD_SCM_CSV:
            sprintf(ConfigPath, "%s%05d.csv", md->SaveDirectory, md->_FileIndex++);
            md->ConfigFilep = f_open(ConfigPath, "w");
            for (int i=0; i < md->potsys.atomkinds_num; i++)
            {
                for (int k=0; k < md->TotalNoOfParticles; k++)
                    if (globaldata[k].atomkind == i)
                        fprintf(md->ConfigFilep, "%.2f, %.2f, %.2f, %d, %.4f\n",
                            globaldata[k].x[0], globaldata[k].x[1],
                            globaldata[k].x[2], i, globaldata[k].var);
            }
            fclose(md->ConfigFilep);
            break;

        case FMD_SCM_VTF:
            save_VTF_file(md, globaldata);
            break;
    }

    free(globaldata);
}

void fmd_matt_setDesiredTemperature(fmd_t *md, fmd_real_t DesiredTemperature)
{
    md->DesiredTemperature = DesiredTemperature;
}
