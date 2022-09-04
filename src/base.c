/*
  base.c: This file is part of Free Molecular Dynamics

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

#include <stdarg.h>
#include <stddef.h>
#include "base.h"
#include "md_ghost.h"
#include "forces.h"
#include "timer.h"
#include "molecule.h"
#include "array.h"
#include "turi.h"
#include "general.h"
#include "cell.h"

typedef struct
{
    float x[DIM];
    float var;
    unsigned atomkind;
} config_atom_t;

typedef struct
{
    fmd_rtuple_t x;
    fmd_rtuple_t v;
    int GroupID;
    unsigned atomkind;
} state_atom_t;

static char formatstr_3xpoint16e[] = "%.16e\t%.16e\t%.16e\n";

void _fmd_initialize_grid(cell_t ***grid, cellinfo_t *cinfo, unsigned dim1, unsigned dim2, unsigned dim3)
{
    for (unsigned i=0; i < dim1; i++)
        for (unsigned j=0; j < dim2; j++)
            for (unsigned k=0; k < dim3; k++)
                _fmd_cell_init(cinfo, &grid[i][j][k]);
}

void compLocOrdParam(fmd_t *md)
{
/*
    float latticeParameter = md->EAM.elements[0].latticeParameter;
    float rCutSqd = sqrr(1.32 * latticeParameter);
    fmd_ituple_t ic, jc, kc;
    int d;
    ParticleListItem_t *item1_p, *item2_p;
    float arg, fmd_real_t, img;
    int Z;
    float r2, rv[3];
    float q[6][3] = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{0,1,1},{1,0,1}};
    int i;

    if (md->LOP_iteration == 0)
        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (item1_p = ARRAY_ELEMENT(md->SubDomain.grid, ic); item1_p != NULL; item1_p = item1_p->next_p)
            {
                item1_p->P.LocOrdParam = 0.;
                for (d=0; d<3; d++)
                    item1_p->P.x_avgd[d] = 0.0;
            }

    (md->LOP_iteration)++;

    for (i=0; i<6; i++)
        for (d=0; d<3; d++)
            q[i][d] *= 4.0 * M_PI / latticeParameter * q[i][d];

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item1_p = ARRAY_ELEMENT(md->SubDomain.grid, ic); item1_p != NULL; item1_p = item1_p->next_p)
        {
            fmd_real_t = img = 0.;
            Z = 0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic[0]-1; kc[0]<=ic[0]+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic[1]-1; kc[1]<=ic[1]+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic[2]-1; kc[2]<=ic[2]+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = ARRAY_ELEMENT(md->SubDomain.grid, jc); item2_p != NULL; item2_p = item2_p->next_p)
                            if (item1_p != item2_p)
                            {
                                for (d=0; d<3; d++)
                                {
                                    if (md->ns[d] == 1)
                                    {
                                        if (kc[d]==-1)
                                            rv[d] = item1_p->P.x[d] - item2_p->P.x[d] + md->l[d];
                                        else
                                            if (kc[d] == md->nc[d])
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d] - md->l[d];
                                            else
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                    }
                                    else
                                        rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                }
                                r2 = sqrr(rv[0])+sqrr(rv[1])+sqrr(rv[2]);
                                if (r2 < rCutSqd)
                                {
                                    Z++;
                                    for (i=0; i<6; i++)
                                    {
                                        arg = 0.;
                                        for (d=0; d<3; d++)
                                            arg += q[i][d] * rv[d];
                                        fmd_real_t += cosf(arg);
                                        img  += sinf(arg);
                                    }
                                }
                            }
                    }
                }
            }
            item1_p->P.LocOrdParam += (sqrr(fmd_real_t)+sqrr(img)) / sqrr(Z);
            for (d=0; d<3; d++)
                item1_p->P.x_avgd[d] += item1_p->P.x[d];
        }

    if (md->LOP_iteration != md->LOP_period) return;

    _fmd_ghostparticles_update_LocOrdParam(md);

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item1_p = ARRAY_ELEMENT(md->SubDomain.grid, ic); item1_p != NULL; item1_p = item1_p->next_p)
        {
            fmd_real_t = item1_p->P.LocOrdParam;
            Z = 0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic[0]-1; kc[0]<=ic[0]+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic[1]-1; kc[1]<=ic[1]+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic[2]-1; kc[2]<=ic[2]+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = ARRAY_ELEMENT(md->SubDomain.grid, jc); item2_p != NULL; item2_p = item2_p->next_p)
                            if (item1_p != item2_p)
                            {
                                for (d=0; d<3; d++)
                                {
                                    if (md->ns[d] == 1)
                                    {
                                        if (kc[d]==-1)
                                            rv[d] = item1_p->P.x[d] - item2_p->P.x[d] + md->l[d];
                                        else
                                            if (kc[d] == md->nc[d])
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d] - md->l[d];
                                            else
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                    }
                                    else
                                        rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                }
                                r2 = sqrr(rv[0])+sqrr(rv[1])+sqrr(rv[2]);
                                if (r2 < rCutSqd)
                                {
                                    Z++;
                                    fmd_real_t += item2_p->P.LocOrdParam;

                                }
                            }
                    }
                }
            }
            item1_p->LocOrdParamAvg = fmd_real_t / ((Z+1) * 36 * md->LOP_period);
            for (d=0; d<3; d++)
                item1_p->P.x_avgd[d] /= md->LOP_period;
        }

    md->LOP_iteration = 0;
*/
}

/* not correct under periodic boundary conditions
   see [J. Chem. Phys. 131, 154107 (2009)] */
static fmd_real_t compVirial_internal(fmd_t *md)
{
    fmd_ituple_t ic;
    fmd_real_t virial = 0.0;
    fmd_real_t virial_global;
    unsigned i;
    cell_t *c;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->SubDomain.grid, ic), i=0; i < c->parts_num; i++)
            virial += POS(c, i, 0) * FRC(c, i, 0) +
                      POS(c, i, 1) * FRC(c, i, 1) +
                      POS(c, i, 2) * FRC(c, i, 2);

    MPI_Reduce(&virial, &virial_global, 1, FMD_MPI_REAL, MPI_SUM, RANK0, md->MD_comm);

    return virial_global;
}

void createCommunicators(fmd_t *md)
{
    int mdnum, i;
    MPI_Group world_group, MD_group;
    int *ranks;

    /* create MD_comm */

    mdnum = md->ns[0] * md->ns[1] * md->ns[2];

    ranks = (int *)m_alloc(mdnum * sizeof(int));

    for (i=0; i<mdnum; i++)
        ranks[i] = i;

    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, mdnum, ranks, &MD_group);
    MPI_Comm_create(MPI_COMM_WORLD, MD_group, &md->MD_comm);
    MPI_Group_free(&world_group);
    MPI_Group_free(&MD_group);

    free(ranks);
}

void findLimits(fmd_t *md, fmd_rtuple_t LowerLimit, fmd_rtuple_t UpperLimit)
{
    fmd_ituple_t ic;
    int d;
    fmd_rtuple_t LocalLower, LocalUpper;

    LocalLower[0] = LocalLower[1] = LocalLower[2] = DBL_MAX;
    LocalUpper[0] = LocalUpper[1] = LocalUpper[2] = DBL_MIN;

    int i;
    cell_t *c;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->SubDomain.grid, ic), i=0; i < c->parts_num; i++)
            for (d=0; d<3; d++)
            {
                if (POS(c, i, d) < LocalLower[d])
                    LocalLower[d] = POS(c, i, d);
                if (POS(c, i, d) > LocalUpper[d])
                    LocalUpper[d] = POS(c, i, d);
            }

    MPI_Allreduce(LocalLower, LowerLimit, 3, FMD_MPI_REAL, MPI_MIN, md->MD_comm);
    MPI_Allreduce(LocalUpper, UpperLimit, 3, FMD_MPI_REAL, MPI_MAX, md->MD_comm);
}

static void identifyProcess(fmd_t *md)
{
    int mdnum;

    mdnum = md->ns[0] * md->ns[1] * md->ns[2];
    if (md->world_rank < mdnum)
        md->Is_MD_process = FMD_TRUE;
    else
        md->Is_MD_process = FMD_FALSE;
}

void fmd_matt_addVelocity(fmd_t *md, int GroupID, fmd_real_t vx, fmd_real_t vy, fmd_real_t vz)
{
    cell_t ***grid;
    const int *start, *stop;
    fmd_ituple_t ic;

    if (md->ParticlesDistributed)
    {
        grid = md->SubDomain.grid;
        start = md->SubDomain.ic_start;
        stop = md->SubDomain.ic_stop;
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
            if (GroupID == -1 || GroupID == c->GroupID[i])
            {
                VEL(c, i, 0) += vx;
                VEL(c, i, 1) += vy;
                VEL(c, i, 2) += vz;
            }
}

/* calculate GroupTemperature, ActiveGroupParticlesNum and TotalKineticEnergy from "global grid" */
static void calculate_GroupTemperature_etc(fmd_t *md)
{
    if (md->Is_MD_comm_root)
    {
        fmd_ituple_t ic;
        fmd_real_t m_vSqd_Sum = 0.0;
        cell_t *c;
        int i;

        md->ActiveGroupParticlesNum = 0;

        LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
            for (c=&ARRAY_ELEMENT(md->global_grid, ic), i=0; i < c->parts_num; i++)
            {
                if (md->ActiveGroup != ACTIVE_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                    continue;

                md->ActiveGroupParticlesNum++;

                fmd_real_t mass = md->potsys.atomkinds[c->atomkind[i]].mass;
                for (int d=0; d<DIM; d++)
                    m_vSqd_Sum += mass * sqrr(VEL(c, i, d));
            }

        md->GroupTemperature = m_vSqd_Sum / (3.0 * md->ActiveGroupParticlesNum * K_BOLTZMANN);
    }

    MPI_Bcast(&md->ActiveGroupParticlesNum, 1, MPI_UNSIGNED, RANK0, md->MD_comm);
    MPI_Bcast(&md->GroupTemperature, 1, FMD_MPI_REAL, RANK0, md->MD_comm);

    md->TotalKineticEnergy = 3.0/2.0 * md->ActiveGroupParticlesNum * K_BOLTZMANN * md->GroupTemperature;
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

    LOOP3D(icl, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {
        _fmd_conv_ic_loc_to_glob(md, icl, icg);

        cl = &ARRAY_ELEMENT(md->SubDomain.grid, icl);
        cg = &ARRAY_ELEMENT(md->global_grid, icg);

        cl->parts_num = cg->parts_num;
        md->SubDomain.NumberOfParticles += cl->parts_num;

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

void unpack_for_matt_distribute(fmd_t *md, void *packbuf, int bufsize)
{
    fmd_ituple_t ic;
    cell_t *c;
    int pos = 0;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {
        c = &ARRAY_ELEMENT(md->SubDomain.grid, ic);

        MPI_Unpack(packbuf, bufsize, &pos, &c->parts_num, 1, MPI_UNSIGNED, md->MD_comm);

        if (c->parts_num > 0)
        {
            md->SubDomain.NumberOfParticles += c->parts_num;

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
    if (md->SubDomain.grid == NULL) fmd_subd_init(md);

    calculate_GroupTemperature_etc(md);

    if (md->Is_MD_comm_root)
    {
        for (int i=1; i < md->SubDomain.numprocs; i++)  /* for all processes execpt the root */
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

    md->GlobalGridExists = FMD_FALSE;
    md->ParticlesDistributed = FMD_TRUE;
}

void fmd_io_loadState(fmd_t *md, fmd_string_t file, fmd_bool_t UseTime)
{
    FILE *fp;
    char name[3];
    fmd_ituple_t ic;
    fmd_real_t StateFileTime;
    unsigned ParticlesNum;
    fmd_real_t l0, l1, l2;
    int PBC0, PBC1, PBC2;

    if (md->Is_MD_comm_root)
    {
        fp = f_open(file, "r");

        assert( fscanf(fp, "%lf", &StateFileTime) == 1 ); /* TO-DO: handle error */

        if (UseTime) md->time = StateFileTime;

        assert( fscanf(fp, "%u\n", &ParticlesNum) == 1 ); /* TO-DO: handle error */
        assert( fscanf(fp, "%lf%lf%lf", &l0, &l1, &l2) == 3 ); /* TO-DO: handle error */
        assert( fscanf(fp, "%d%d%d", &PBC0, &PBC1, &PBC2) == 3 ); /* TO-DO: handle error */
    }

    if (!md->BoxSizeDetermined)
    {
        if (md->Is_MD_comm_root)
        {
            md->l[0] = l0;
            md->l[1] = l1;
            md->l[2] = l2;
        }

        MPI_Bcast(&md->l, 3, FMD_MPI_REAL, RANK0, md->MD_comm);
        md->BoxSizeDetermined = FMD_TRUE;
    }

    if (!md->PBCdetermined)
    {
        if (md->Is_MD_comm_root)
        {
            md->PBC[0] = PBC0;
            md->PBC[1] = PBC1;
            md->PBC[2] = PBC2;
        }

        MPI_Bcast(&md->PBC, 3, MPI_INT, RANK0, md->MD_comm);
        md->PBCdetermined = FMD_TRUE;
    }

    if (!md->GlobalGridExists)
        fmd_box_createGrid(md, md->CutoffRadius);

    if (UseTime)
        MPI_Bcast(&md->time, 1, FMD_MPI_REAL, RANK0, md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        for (unsigned i=0; i < ParticlesNum; i++)
        {
            int GroupID;
            unsigned atomkind;

            assert( fscanf(fp, "%s%d", name, &GroupID) == 2 ); /* TO-DO: handle error */
            assert(GroupID >= 0);  /* TO-DO: handle error */

            unsigned j;

            for (j=0; j < md->potsys.atomkinds_num; j++)
                if (strcmp(name, md->potsys.atomkinds[j].name) == 0)
                {
                    atomkind = j;
                    break;
                }

            /* TO-DO: what if the name doesn't exist in potsys? */
            assert(j < md->potsys.atomkinds_num);

            fmd_rtuple_t x, v;

            assert( fscanf(fp, "%lf%lf%lf", &x[0], &x[1], &x[2]) == 3 ); /* TO-DO: handle error */
            assert( fscanf(fp, "%lf%lf%lf", &v[0], &v[1], &v[2]) == 3 ); /* TO-DO: handle error */

            for (int d=0; d<DIM; d++)
                ic[d] = (int)floor(x[d] / md->cellh[d]);

            cell_t *c;

            c = &ARRAY_ELEMENT(md->global_grid, ic);

            j = _fmd_cell_new_particle(md, c);

            c->GroupID[j] = GroupID;
            c->atomkind[j] = atomkind;
            c->AtomID[j] = md->TotalNoOfParticles++;

            for (int d=0; d<DIM; d++)
            {
                POS(c, j, d) = x[d];
                VEL(c, j, d) = v[d];
            }
        }

        fclose(fp);
    }
}

void _fmd_refreshGrid(fmd_t *md)
{
    fmd_ituple_t ic, jc;

    /* iterate over all cells */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {
        /* iterate over all particles in cell ic */

        cell_t *c = &ARRAY_ELEMENT(md->SubDomain.grid, ic);
        int i = 0;  /* particle index */

        while (i < c->parts_num)
        {
            for (int d=0; d<DIM; d++)
            {
                jc[d] = (int)floor(POS(c, i, d) / md->cellh[d]) - md->SubDomain.ic_global_firstcell[d]
                        + md->SubDomain.ic_start[d];

                if (jc[d] < 0 || jc[d] >= md->SubDomain.cell_num[d])
                {
                    fprintf(stderr, "ERROR: Unexpected particle position!\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR_UNEXPECTED_PARTICLE_POSITION);
                }
            }

            if ((ic[0] != jc[0]) || (ic[1] != jc[1]) || (ic[2] != jc[2]))
            {
                cell_t *c2 = &ARRAY_ELEMENT(md->SubDomain.grid, jc);
                unsigned j = _fmd_cell_new_particle(md, c2);

                _fmd_cell_copy_atom_from_cell_to_cell(c, i, c2, j);
                _fmd_cell_remove_atom(md, c, i);
            }
            else
                i++;
        }
    }

    /* now particles in ghost cells migrate to neighbour subdomains */
    _fmd_particles_migrate(md);
}

void rescaleVelocities(fmd_t *md)
{
    fmd_ituple_t ic;

    fmd_real_t scale = sqrt(md->DesiredTemperature / md->GroupTemperature);

    cell_t *cell;
    int i;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &ARRAY_ELEMENT(md->SubDomain.grid, ic), i=0; i < cell->parts_num; i++)
            for (int d=0; d<3; d++)
                VEL(cell, i, d) *= scale;

    md->GroupTemperature = md->DesiredTemperature;
}

static config_atom_t *prepare_localdata_for_saveconfig(fmd_t *md)
{
    config_atom_t *localdata;
    fmd_ituple_t ic;
    cell_t *c;
    int pind;
    int k = 0;

    localdata = (config_atom_t *)m_alloc(md->SubDomain.NumberOfParticles * sizeof(config_atom_t));

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->SubDomain.grid, ic), pind=0; pind < c->parts_num; pind++)
        {
            if (md->CompLocOrdParam)
            {
//                 localdata[k].x[0] = (float)pc->x_avgd[0];
//                 localdata[k].x[1] = (float)pc->x_avgd[1];
//                 localdata[k].x[2] = (float)pc->x_avgd[2];
            }
            else
            {
                localdata[k].x[0] = (float)POS(c, pind, 0);
                localdata[k].x[1] = (float)POS(c, pind, 1);
                localdata[k].x[2] = (float)POS(c, pind, 2);
            }
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

    if (md->Is_MD_comm_root) nums = (unsigned *)m_alloc(md->SubDomain.numprocs * sizeof(unsigned));
    MPI_Gather(&md->SubDomain.NumberOfParticles, 1, MPI_UNSIGNED, nums, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        md->TotalNoOfParticles = 0;

        for (int i=0; i < md->SubDomain.numprocs; i++)
            md->TotalNoOfParticles += nums[i];

        int displ = 0;

        globaldata = (config_atom_t *)m_alloc(md->TotalNoOfParticles * sizeof(config_atom_t));
        displs     = (int *)m_alloc(md->SubDomain.numprocs * sizeof(int));

        for (int k=0; k < md->SubDomain.numprocs; k++)
        {
            displs[k] = displ;
            displ += nums[k];
        }
    }

    MPI_Gatherv(localdata, md->SubDomain.NumberOfParticles, md->mpi_types.mpi_configa,
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
        case FMD_SCM_XYZ_PARTICLESNUM:
            if (md->TotalNoOfParticles != md->_OldNumberOfParticles)
            {
                if (md->_OldNumberOfParticles != -1) fclose(md->ConfigFilep);
                sprintf(ConfigPath, "%s%d.xyz", md->SaveDirectory, md->TotalNoOfParticles);
                md->ConfigFilep = f_open(ConfigPath, "w");
                md->_OldNumberOfParticles = md->TotalNoOfParticles;
            }
            save_XYZ_data(md, globaldata);
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

void fmd_io_saveState(fmd_t *md, fmd_string_t filename)
{
    unsigned *nums;
    fmd_ituple_t ic;
    char StateFilePath[MAX_PATH_LENGTH];
    FILE *fp;
    MPI_Status status;
    cell_t *c;
    unsigned pi;

    if (md->Is_MD_comm_root)
    {
        nums = (unsigned *)m_alloc(md->SubDomain.numprocs * sizeof(unsigned));

        MPI_Gather(&md->SubDomain.NumberOfParticles, 1, MPI_UNSIGNED, nums, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

        md->TotalNoOfParticles = 0;

        for (int k=0; k < md->SubDomain.numprocs; k++)
            md->TotalNoOfParticles += nums[k];

        sprintf(StateFilePath, "%s%s", md->SaveDirectory, filename);
        fp = f_open(StateFilePath, "w");

        fprintf(fp, "%.16e\n", md->time);
        fprintf(fp, "%u\n", md->TotalNoOfParticles);
        fprintf(fp, formatstr_3xpoint16e, md->l[0], md->l[1], md->l[2]);
        fprintf(fp, "%d %d %d\n", md->PBC[0], md->PBC[1], md->PBC[2]);

        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (c = &ARRAY_ELEMENT(md->SubDomain.grid, ic), pi=0; pi < c->parts_num; pi++)
            {
                fprintf(fp, "%s %d\n", md->potsys.atomkinds[c->atomkind[pi]].name, c->GroupID[pi]);
                fprintf(fp, formatstr_3xpoint16e, POS(c, pi, 0), POS(c, pi, 1), POS(c, pi, 2));
                fprintf(fp, formatstr_3xpoint16e, VEL(c, pi, 0), VEL(c, pi, 1), VEL(c, pi, 2));
            }

        for (int i=1; i < md->SubDomain.numprocs; i++)
        {
            state_atom_t *states = (state_atom_t *)m_alloc(nums[i] * sizeof(state_atom_t));

            MPI_Recv(states, nums[i], md->mpi_types.mpi_statea, i, 150, md->MD_comm, &status);

            for (unsigned k=0; k < nums[i]; k++)
            {
                state_atom_t *s = states + k;

                fprintf(fp, "%s %d\n", md->potsys.atomkinds[s->atomkind].name, s->GroupID);
                fprintf(fp, formatstr_3xpoint16e, s->x[0], s->x[1], s->x[2]);
                fprintf(fp, formatstr_3xpoint16e, s->v[0], s->v[1], s->v[2]);
            }

            free(states);
        }

        free(nums);
        fclose(fp);
    }
    else
    {
        MPI_Gather(&md->SubDomain.NumberOfParticles, 1, MPI_UNSIGNED, nums, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

        state_atom_t *states = (state_atom_t *)m_alloc(md->SubDomain.NumberOfParticles * sizeof(state_atom_t));

        unsigned k = 0;

        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (c = &ARRAY_ELEMENT(md->SubDomain.grid, ic), pi=0; pi < c->parts_num; pi++)
            {
                for (int d=0; d<DIM; d++)
                    states[k].x[d] = POS(c, pi, d);

                for (int d=0; d<DIM; d++)
                    states[k].v[d] = VEL(c, pi, d);

                states[k].GroupID = c->GroupID[pi];
                states[k].atomkind = c->atomkind[pi];

                k++;
            }

        MPI_Send(states, md->SubDomain.NumberOfParticles, md->mpi_types.mpi_statea, RANK0, 150, md->MD_comm);

        free(states);
    }
}

void fmd_matt_setDesiredTemperature(fmd_t *md, fmd_real_t DesiredTemperature)
{
    md->DesiredTemperature = DesiredTemperature;
}

void fmd_box_setPBC(fmd_t *md, fmd_bool_t PBCx, fmd_bool_t PBCy, fmd_bool_t PBCz)
{
    md->PBC[0] = PBCx;
    md->PBC[1] = PBCy;
    md->PBC[2] = PBCz;
    md->PBCdetermined = FMD_TRUE;
}

void fmd_box_setSubDomains(fmd_t *md, int dimx, int dimy, int dimz)
{
    md->ns[0] = dimx;
    md->ns[1] = dimy;
    md->ns[2] = dimz;
    identifyProcess(md);
    createCommunicators(md);
    if (md->Is_MD_process)
    {
        MPI_Comm_size(md->MD_comm, &md->SubDomain.numprocs);
        MPI_Comm_rank(md->MD_comm, &md->SubDomain.myrank);
        if (md->SubDomain.myrank == RANK0)
            md->Is_MD_comm_root = FMD_TRUE;
    }
}

static void create_mpi_types(fmd_t *md)
{
    MPI_Datatype temptype;

    /* create and commit mpi_ituple */
    MPI_Type_contiguous(DIM, MPI_INT, &md->mpi_types.mpi_ituple);
    MPI_Type_commit(&md->mpi_types.mpi_ituple);

    /* create and commit mpi_rtuple */
    MPI_Type_contiguous(DIM, FMD_MPI_REAL, &md->mpi_types.mpi_rtuple);
    MPI_Type_commit(&md->mpi_types.mpi_rtuple);

    /* create and commit mpi_configa */
    int configblocklen[3] = {DIM, 1, 1};
    MPI_Aint configdisplc[3] = {offsetof(config_atom_t, x),
                                offsetof(config_atom_t, var),
                                offsetof(config_atom_t, atomkind)};
    MPI_Datatype configtype[3] = {MPI_FLOAT, MPI_FLOAT, MPI_UNSIGNED};
    MPI_Type_create_struct(3, configblocklen, configdisplc, configtype, &temptype);
    MPI_Type_create_resized(temptype, 0, sizeof(config_atom_t), &md->mpi_types.mpi_configa);
    MPI_Type_free(&temptype);
    MPI_Type_commit(&md->mpi_types.mpi_configa);

    /* create and commit mpi_statea */
    int stateblocklen[4] = {1, 1, 1, 1};
    MPI_Aint statedisplc[4] = {offsetof(state_atom_t, x),
                               offsetof(state_atom_t, v),
                               offsetof(state_atom_t, GroupID),
                               offsetof(state_atom_t, atomkind)};
    MPI_Datatype statetype[4] = {md->mpi_types.mpi_rtuple,
                                 md->mpi_types.mpi_rtuple,
                                 MPI_INT,
                                 MPI_UNSIGNED};
    MPI_Type_create_struct(4, stateblocklen, statedisplc, statetype, &temptype);
    MPI_Type_create_resized(temptype, 0, sizeof(state_atom_t), &md->mpi_types.mpi_statea);
    MPI_Type_free(&temptype);
    MPI_Type_commit(&md->mpi_types.mpi_statea);
}

static void free_mpi_types(fmd_t *md)
{
    MPI_Type_free(&md->mpi_types.mpi_ituple);
    MPI_Type_free(&md->mpi_types.mpi_rtuple);
    MPI_Type_free(&md->mpi_types.mpi_configa);
    MPI_Type_free(&md->mpi_types.mpi_statea);
}

fmd_t *fmd_create()
{
    fmd_t *md = (fmd_t *)m_alloc(sizeof(fmd_t));

    md->MPI_initialized_by_me = FMD_FALSE;

    int isMPIInitialized;

    MPI_Initialized(&isMPIInitialized);

    if (!isMPIInitialized)
    {
        MPI_Init(NULL, NULL);
        md->MPI_initialized_by_me = FMD_TRUE;
    }

    omp_set_num_threads(1);

    MPI_Comm_size(MPI_COMM_WORLD, &(md->world_numprocs));
    MPI_Comm_rank(MPI_COMM_WORLD, &(md->world_rank));
    md->LOP_iteration = 0;
    md->UseAutoStep = FMD_FALSE;
    md->time = 0.0;
    md->time_iteration = 0;
    md->SaveDirectory[0] = '\0';
    md->CompLocOrdParam = FMD_FALSE;
    md->SubDomain.grid = NULL;
    md->TotalNoOfParticles = 0;
    md->TotalNoOfMolecules = 0;
    md->ActiveGroup = ACTIVE_GROUP_ALL;             /* all groups are active by default */
    md->ParticlesDistributed = FMD_FALSE;
    md->GlobalGridExists = FMD_FALSE;
    md->global_grid = NULL;
    md->BoxSizeDetermined = FMD_FALSE;
    md->PBCdetermined = FMD_FALSE;
    md->Is_MD_comm_root = FMD_FALSE;
    md->EventHandler = NULL;
    md->timers = NULL;
    md->timers_num = 0;
    md->turies = NULL;
    md->turies_num = 0;
    md->active_ttm_turi = NULL;
    md->DesiredTemperature = 300.0;
    md->SaveConfigMode = FMD_SCM_XYZ_PARTICLESNUM;
    md->cell_increment = 10;
    md->_OldNumberOfParticles = -1;
    md->_FileIndex = 0;
    md->_OldTotalMDEnergy = 0.0;
    md->_PrevFailedMDEnergy = 0.0;
    fmd_potsys_init(md);
    create_mpi_types(md);
    _fmd_h5_ds_init(&md->h5_dataspaces);

    // this must be the last statement before return
    md->WallTimeOrigin = MPI_Wtime();

    return md;
}

void fmd_box_setSize(fmd_t *md, fmd_real_t sx, fmd_real_t sy, fmd_real_t sz)
{
    if (!md->GlobalGridExists)
    {
        md->l[0] = sx;
        md->l[1] = sy;
        md->l[2] = sz;
        md->BoxSizeDetermined = FMD_TRUE;
    }
}

fmd_real_t fmd_proc_getWallTime(fmd_t *md)
{
    return (MPI_Wtime() - md->WallTimeOrigin);
}

fmd_bool_t fmd_proc_isMD(fmd_t *md)
{
    return md->Is_MD_process;
}

fmd_bool_t fmd_proc_isRoot(fmd_t *md)
{
    return md->Is_MD_comm_root;
}

void fmd_box_createGrid(fmd_t *md, fmd_real_t cutoff)
{
    for (int d=0; d<3; d++)
    {
        md->nc[d] = (int)(md->l[d] / cutoff);
        if ((md->nc[d] < 3) && md->PBC[d])
        {
            fprintf(stderr, "ERROR: nc[%d] = %d. Under PBC, this must be greater than 2!\n", d, md->nc[d]);
            MPI_Abort(MPI_COMM_WORLD, ERROR_NCELL_TOO_SMALL);
        }
        md->cellh[d] = md->l[d] / md->nc[d];
    }

    _fmd_cellinfo_init(&md->cellinfo);

    if (md->Is_MD_comm_root)
    {
        md->global_grid = (cell_t ***)_fmd_array_ordinary3d_create(md->nc, sizeof(cell_t));
        assert(md->global_grid != NULL);
        /* TO-DO: handle memory error */

        _fmd_initialize_grid(md->global_grid, &md->cellinfo, md->nc[0], md->nc[1], md->nc[2]);
    }

    md->GlobalGridExists = FMD_TRUE;
    md->CutoffRadius = cutoff;
}

void fmd_io_setSaveDirectory(fmd_t *md, fmd_string_t directory)
{
    strcpy(md->SaveDirectory, directory);
}

void fmd_io_setSaveConfigMode(fmd_t *md, fmd_SaveConfigMode_t mode)
{
    md->SaveConfigMode = mode;
}

void fmd_io_printf(fmd_t *md, const fmd_string_t restrict format, ...)
{
    if (md->Is_MD_process && md->Is_MD_comm_root)
    {
        va_list argptr;

        va_start(argptr, format);
        vprintf(format, argptr);
        va_end(argptr);
    }
}

fmd_real_t fmd_matt_getTotalEnergy(fmd_t *md)
{
    return md->TotalKineticEnergy + md->TotalPotentialEnergy;
}

void fmd_matt_giveTemperature(fmd_t *md, int GroupID)
{
    cell_t ***grid;
    const int *start, *stop;
    fmd_ituple_t ic;

    if (md->ParticlesDistributed)
    {
        grid = md->SubDomain.grid;
        start = md->SubDomain.ic_start;
        stop = md->SubDomain.ic_stop;
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
            if (GroupID == -1 || GroupID == c->GroupID[pi])
            {
                fmd_real_t mass = md->potsys.atomkinds[c->atomkind[pi]].mass;
                fmd_real_t StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

                for (int d=0; d<DIM; d++)
                    VEL(c, pi, d) = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
            }

    gsl_rng_free(rng);
}

fmd_real_t fmd_matt_getGroupTemperature(fmd_t *md)
{
    return md->GroupTemperature;
}

void fmd_free(fmd_t *md)
{
    fmd_subd_free(md);
    fmd_potsys_free(md);
    fmd_timer_free(md);
    free_mpi_types(md);
    _fmd_h5_ds_free(&md->h5_dataspaces);
    fmd_turi_free(md);
    free(md);

    if (md->MPI_initialized_by_me)
    {
        int Is_MPI_Finalized;
        MPI_Finalized(&Is_MPI_Finalized);
        if (!Is_MPI_Finalized) MPI_Finalize();
    }
}

void fmd_setEventHandler(fmd_t *md, fmd_EventHandler_t func)
{
    md->EventHandler = func;
}
