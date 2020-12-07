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
#include "base.h"
#include "md_ghost.h"
#include "forces.h"
#ifdef USE_TTM
#include "ttm.h"
#endif
#include "timer.h"
#include "molecule.h"
#include "array.h"
#include "turi.h"
#include "general.h"

static void refreshGrid(fmd_t *md, int reverse);

void _fmd_initialize_grid(cell_t ***grid, unsigned dim1, unsigned dim2, unsigned dim3)
{
     for (int i=0; i < dim1; i++)
        for (int j=0; j < dim2; j++)
            for (int k=0; k < dim3; k++)
            {
                grid[i][j][k].capacity = 0;
                grid[i][j][k].parts = NULL;
                grid[i][j][k].parts_num = 0;
            }
}

void _fmd_cleanGridSegment(cell_t ***grid, fmd_ituple_t ic_from, fmd_ituple_t ic_to)
{
    fmd_ituple_t ic;

    LOOP3D(ic, ic_from, ic_to)
        FREE_CELL(grid[ic[0]][ic[1]][ic[2]]);
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
            for (item1_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
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
        for (item1_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
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
                        for (item2_p = md->SubDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
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

    fmd_ghostparticles_update_LocOrdParam(md);

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item1_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
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
                        for (item2_p = md->SubDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
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

void fmd_dync_VelocityVerlet_startStep(fmd_t *md, fmd_bool_t UseThermostat)
{
    fmd_ituple_t ic;
    int d;
    fmd_real_t velocityScale, mass;
    fmd_real_t x;

    if (UseThermostat) velocityScale = sqrt(1 + md->delta_t / md->BerendsenThermostatParam *
                       (md->DesiredTemperature / md->GlobalTemperature - 1));

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {
        /* iterate over all particles in cell ic */

        int i=0;
        cell_t *cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]];

        while (i < cell->parts_num)
        {
            particle_t *p = &cell->parts[i];

            if (!(md->ActiveGroup == -1 || p->core.GroupID == md->ActiveGroup))
            {
                i++;
                continue;
            }

            mass = md->potsys.atomkinds[p->core.atomkind].mass;

            for (d=0; d<3; d++)
            {
                /*if (md->UseAutoStep)
                {
                    p->core.v_bak[d] = p->core.v[d];
                    p->core.x_bak[d] = p->core.x[d];
                }*/

                if (UseThermostat) p->core.v[d] *= velocityScale;

                p->core.v[d] += md->delta_t * 0.5 / mass * p->F[d];
                x = p->core.x[d] + md->delta_t * p->core.v[d];

                if ( (md->ns[d] == 1) && ((x < 0.0) || (x >= md->l[d])) )
                {
                    if (!md->PBC[d])
                    {
                        REMOVE_PARTICLE_FROM_CELL(*cell, i);
                        (md->SubDomain.NumberOfParticles)--;
                        i--;
                        break;
                    }
                    else
                        if (x < 0.0) x += md->l[d]; else x -= md->l[d];
                }

                p->core.x[d] = x;
            }

            i++;
        }
    }

    refreshGrid(md, 0);
}

int fmd_dync_VelocityVerlet_finishStep(fmd_t *md)
{
    fmd_ituple_t ic;
    int d;
    fmd_real_t m_vSqd_Sum = 0, m_vSqd_SumSum;
    fmd_real_t mass;
    int returnVal = 0;
    int particlesNum = 0;
    fmd_rtuple_t momentumSum = {0., 0., 0.};
    cell_t *cell;
    int i;

    for (ic[0] = md->SubDomain.ic_start[0]; ic[0] < md->SubDomain.ic_stop[0]; ic[0]++)
    {
        for (ic[1] = md->SubDomain.ic_start[1]; ic[1] < md->SubDomain.ic_stop[1]; ic[1]++)
            for (ic[2] = md->SubDomain.ic_start[2]; ic[2] < md->SubDomain.ic_stop[2]; ic[2]++)
                for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
                {
                    particle_t *p = &cell->parts[i];

                    if (!(md->ActiveGroup == -1 || p->core.GroupID == md->ActiveGroup))
                        continue;
                    particlesNum++;
                    mass = md->potsys.atomkinds[p->core.atomkind].mass;
                    for (d=0; d<3; d++)
                    {
                        p->core.v[d] += md->delta_t * 0.5 / mass * p->F[d];
                        momentumSum[d] += mass * p->core.v[d];
                    }
                    m_vSqd_Sum += mass * ( sqrr(p->core.v[0]) +
                                           sqrr(p->core.v[1]) +
                                           sqrr(p->core.v[2]) );
                }
    }
    MPI_Reduce(momentumSum, md->TotalMomentum, 3, FMD_MPI_REAL, MPI_SUM,
               ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

    MPI_Allreduce(&particlesNum, &(md->ActiveGroupParticlesNum), 1, MPI_INT, MPI_SUM, md->MD_comm);

    if (md->ActiveGroup == -1)
        md->TotalNoOfParticles = md->ActiveGroupParticlesNum;

    MPI_Allreduce(&m_vSqd_Sum, &m_vSqd_SumSum, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
    md->TotalKineticEnergy = 0.5 * m_vSqd_SumSum;
    md->TotalMDEnergy = md->TotalKineticEnergy + md->TotalPotentialEnergy;

    /*if (md->UseAutoStep)
    {
        if (md->_OldTotalMDEnergy!=0. && fabs((md->TotalMDEnergy-md->_OldTotalMDEnergy)/md->_OldTotalMDEnergy) > md->AutoStepSensitivity)
        {
            if (md->_PrevFailedMDEnergy!=0. && fabs((md->TotalMDEnergy-md->_PrevFailedMDEnergy)/md->_PrevFailedMDEnergy) < md->AutoStepSensitivity)
            {  // was jump in energy due to escape of some energetic particle(s) from simulation box?
                if (md->Is_MD_comm_root)
                    printf("Maybe the jump was caused by departure of some energetic particle(s). Increasing time step...\n");
                returnVal = 2;
            }
            else
            {
                if (md->Is_MD_comm_root)
                {
                    printf("current delta_t = %e\n", md->delta_t);
                    printf("Jump in total MD energy (old=%e new=%e)! Decreasing time step...\n", md->_OldTotalMDEnergy, md->TotalMDEnergy);
                }
                md->_PrevFailedMDEnergy = md->TotalMDEnergy;
                return 1;
            }
        }
        md->_PrevFailedMDEnergy = 0.;
        md->_OldTotalMDEnergy = md->TotalMDEnergy;
    }*/

    md->GlobalTemperature = m_vSqd_SumSum / (3.0 * md->ActiveGroupParticlesNum * K_BOLTZMANN);

    return returnVal;
}

/* not correct under periodic boundary conditions
   see [J. Chem. Phys. 131, 154107 (2009)] */
static fmd_real_t compVirial_internal(fmd_t *md)
{
    fmd_ituple_t ic;
    fmd_real_t virial = 0.0;
    fmd_real_t virial_global;
    int i;
    cell_t *cell;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {
            particle_t *p = &cell->parts[i];

            virial += p->core.x[0] * p->F[0] +
                      p->core.x[1] * p->F[1] +
                      p->core.x[2] * p->F[2];
        }

    MPI_Reduce(&virial, &virial_global, 1, FMD_MPI_REAL, MPI_SUM,
      ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

    return virial_global;
}

void createCommunicators(fmd_t *md)
{
    int mdnum, i;
    MPI_Group world_group, MD_group;
    int *ranks;

    /* create MD_comm */
    mdnum = md->ns[0] * md->ns[1] * md->ns[2];
    ranks = (int *)malloc(mdnum * sizeof(int));
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
    cell_t *cell;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
            for (d=0; d<3; d++)
            {
                particle_core_t *pc = &cell->parts[i].core;

                if (pc->x[d] < LocalLower[d])
                    LocalLower[d] = pc->x[d];
                if (pc->x[d] > LocalUpper[d])
                    LocalUpper[d] = pc->x[d];
            }
    MPI_Allreduce(LocalLower, LowerLimit, 3, FMD_MPI_REAL, MPI_MIN, md->MD_comm);
    MPI_Allreduce(LocalUpper, UpperLimit, 3, FMD_MPI_REAL, MPI_MAX, md->MD_comm);
}

void handleFileOpenError(FILE *fp, char *filename)
{
    if (fp == NULL)
    {
        fprintf(stderr, "ERROR: Unable to open %s!\n", filename);
        MPI_Abort(MPI_COMM_WORLD, ERROR_UNABLE_OPEN_FILE);
    }
}

void identifyProcess(fmd_t *md)
{
    int mdnum;

    mdnum = md->ns[0] * md->ns[1] * md->ns[2];
    if (md->world_rank < mdnum)
        md->Is_MD_process = FMD_TRUE;
    else
        md->Is_MD_process = FMD_FALSE;
#ifdef USE_TTM
    if (ttm_useExtended && md->world_rank == mdnum)
        ttm_is_extended_process = 1;
    else
        ttm_is_extended_process = 0;
#endif
}

void fmd_matt_setActiveGroup(fmd_t *md, int GroupID)
{
    md->ActiveGroup = GroupID;
}

void fmd_matt_addVelocity(fmd_t *md, int GroupID, fmd_real_t vx, fmd_real_t vy, fmd_real_t vz)
{
    cell_t ***grid;
    int *start, *stop;
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
    cell_t *cell;
    particle_core_t *pc;

    LOOP3D(ic, start, stop)
        for (cell = &grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
            if (pc = &cell->parts[i].core, GroupID == -1 || GroupID == pc->GroupID)
            {
                pc->v[0] += vx;
                pc->v[1] += vy;
                pc->v[2] += vz;
            }
}

void fmd_matt_distribute(fmd_t *md)
{
    int i, k, d, nct, sum_length;
    fmd_ituple_t ic;
    int *ic_length;
    particle_core_t *is_partcores;
    cell_t *cell;
    int pind;

    if (md->SubDomain.grid == NULL) fmd_subd_init(md);

    if (md->Is_MD_comm_root)
    {
        int r, w;
        fmd_ituple_t is, global_icstart, global_icstop;
        fmd_real_t m_vSqd_Sum = 0.0;
        fmd_real_t mass;

#ifdef USE_TTM
        ttm_comp_min_atomsNo(global_grid, s_p);
#endif

        LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
            for (cell=&md->global_grid[ic[0]][ic[1]][ic[2]], pind=0; pind < cell->parts_num; pind++)
            {
                particle_core_t *pc = &cell->parts[pind].core;

                mass = md->potsys.atomkinds[pc->atomkind].mass;
                for (d=0; d<3; d++)
                    m_vSqd_Sum += mass * sqrr(pc->v[d]);
            }

        md->GlobalTemperature = m_vSqd_Sum / (3.0 * md->TotalNoOfParticles * K_BOLTZMANN);

        for (i=0; i < ROOTPROCESS(md->SubDomain.numprocs); i++)
        {
            INVERSEINDEX(i, md->ns, is);
            nct = 1;
            for (d=0; d<3; d++)
            {
                r = md->nc[d] % md->ns[d];
                w = md->nc[d] / md->ns[d];
                if (is[d] < r)
                {
                    global_icstart[d] = is[d] * (w + 1);
                    global_icstop[d] = global_icstart[d] + w + 1;
                    nct *= w + 1;
                }
                else
                {
                    global_icstart[d] = is[d] * w + r;
                    global_icstop[d] = global_icstart[d] + w;
                    nct *= w;
                }
            }
            ic_length = (int *)malloc((nct+1) * sizeof(int));
            k = sum_length = 0;
            LOOP3D(ic, global_icstart, global_icstop)
            {
                ic_length[k] = md->global_grid[ic[0]][ic[1]][ic[2]].parts_num;
                sum_length += ic_length[k++];
            }
            ic_length[nct] = sum_length;

            MPI_Send(ic_length, nct+1, MPI_INT, i, 50, md->MD_comm);

            free(ic_length);

            sum_length *= sizeof(particle_core_t);
            is_partcores = (particle_core_t *)malloc(sum_length);

            k = 0;
            LOOP3D(ic, global_icstart, global_icstop)
                for (cell=&md->global_grid[ic[0]][ic[1]][ic[2]], pind=0; pind < cell->parts_num; pind++)
                    is_partcores[k++] = cell->parts[pind].core;

            MPI_Send(is_partcores, sum_length, MPI_CHAR, i, 51, md->MD_comm);

            free(is_partcores);
        }

        md->SubDomain.NumberOfParticles = 0;
        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        {
            cell = &md->global_grid[ ic[0] - md->SubDomain.ic_start[0] + md->SubDomain.ic_global_firstcell[0] ]
                                   [ ic[1] - md->SubDomain.ic_start[1] + md->SubDomain.ic_global_firstcell[1] ]
                                   [ ic[2] - md->SubDomain.ic_start[2] + md->SubDomain.ic_global_firstcell[2] ];
            for (pind = 0; pind < cell->parts_num; pind++)
            {

                cell->parts[pind].neighbors = NULL;
                INSERT_PARTICLE_IN_CELL(cell->parts[pind], md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                ++(md->SubDomain.NumberOfParticles);
            }

            FREE_CELL(*cell);
        }

        _fmd_array_ordinary3d_free((void ***)md->global_grid, md->nc[0], md->nc[1]);
    }
    else
    {
        MPI_Status status;
        int kreceive;

#ifdef USE_TTM
        ttm_comp_min_atomsNo(NULL, s_p);
#endif
        nct = 1;
        for (d=0; d<3; d++)
            nct *= md->SubDomain.cell_num_nonmarg[d];
        ic_length = (int *)malloc((nct+1) * sizeof(int));
        MPI_Recv(ic_length, nct+1, MPI_INT, ROOTPROCESS(md->SubDomain.numprocs),
                 50, md->MD_comm, &status);
        md->SubDomain.NumberOfParticles = sum_length = ic_length[nct];
        sum_length *= sizeof(particle_core_t);
        is_partcores = (particle_core_t *)malloc(sum_length);

        MPI_Recv(is_partcores, sum_length, MPI_CHAR,
                 ROOTPROCESS(md->SubDomain.numprocs), 51, md->MD_comm, &status);

        kreceive = k = 0;
        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        {
            cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]];

            for (i=0; i < ic_length[kreceive]; i++)
            {
                unsigned pind = new_particle(cell);
                cell->parts[pind].core = is_partcores[k++];
                cell->parts[pind].neighbors = NULL;
            }
            kreceive++;
        }

        free(ic_length);
        free(is_partcores);
    }

    MPI_Bcast(&md->TotalNoOfParticles, 1, MPI_UNSIGNED, ROOTPROCESS(md->SubDomain.numprocs),
              md->MD_comm);
    MPI_Bcast(&md->TotalNoOfMolecules, 1, MPI_UNSIGNED, ROOTPROCESS(md->SubDomain.numprocs),
              md->MD_comm);
    MPI_Bcast(&md->GlobalTemperature, 1, FMD_MPI_REAL, ROOTPROCESS(md->SubDomain.numprocs),
              md->MD_comm);

    if (md->TotalNoOfMolecules > 0) _fmd_matt_updateAtomNeighbors(md);

    md->TotalKineticEnergy = 3.0/2.0 * md->TotalNoOfParticles * K_BOLTZMANN * md->GlobalTemperature;
    md->GlobalGridExists = FMD_FALSE;
    md->ParticlesDistributed = FMD_TRUE;
}

void fmd_io_loadState(fmd_t *md, fmd_string_t file, fmd_bool_t UseTime)
{
    FILE *fp;
    char name[3];
    int i, j, d;
    fmd_ituple_t ic;
    fmd_real_t StateFileTime;
    int ParticlesNum;
    fmd_real_t l0, l1, l2;
    int PBC0, PBC1, PBC2;

    if (md->Is_MD_comm_root)
    {
        fp = fopen(file, "r");
        handleFileOpenError(fp, file);
        fscanf(fp, "%lf", &StateFileTime);
        if (UseTime)
            md->MD_time = StateFileTime;
        fscanf(fp, "%d\n", &ParticlesNum);
        md->TotalNoOfParticles += ParticlesNum;
        fscanf(fp, "%lf%lf%lf", &l0, &l1, &l2);
        fscanf(fp, "%d%d%d", &PBC0, &PBC1, &PBC2);
    }

    if (!md->BoxSizeDetermined)
    {
        if (md->Is_MD_comm_root)
        {
            md->l[0] = l0;
            md->l[1] = l1;
            md->l[2] = l2;
        }
        MPI_Bcast(&md->l, 3, FMD_MPI_REAL, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        md->BoxSizeDetermined = 1;
    }

    if (!md->PBCdetermined)
    {
        if (md->Is_MD_comm_root)
        {
            md->PBC[0] = PBC0;
            md->PBC[1] = PBC1;
            md->PBC[2] = PBC2;
        }
        MPI_Bcast(&md->PBC, 3, MPI_INT, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        md->PBCdetermined = 1;
    }

    if (!md->GlobalGridExists)
        fmd_box_createGrid(md, md->CutoffRadius);

    if (UseTime)
        MPI_Bcast(&md->MD_time, 1, FMD_MPI_REAL, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        for (i=0; i < ParticlesNum; i++)
        {
            particle_core_t pc;

            fscanf(fp, "%s%d", name, &pc.GroupID);

            for (j=0; j < md->potsys.atomkinds_num; j++)
                if (strcmp(name, md->potsys.atomkinds[j].name) == 0)
                {
                    pc.atomkind = j;
                    break;
                }

            /* TO-DO: what if the name doesn't exist in potsys? */
            fscanf(fp, "%lf%lf%lf", &pc.x[0], &pc.x[1], &pc.x[2]);
            fscanf(fp, "%lf%lf%lf", &pc.v[0], &pc.v[1], &pc.v[2]);

            for (d=0; d<3; d++)
                ic[d] = (int)floor(pc.x[d] / md->cellh[d]);

            INSERT_PART_CORE_IN_CELL(pc, md->global_grid[ic[0]][ic[1]][ic[2]]);
        }
        fclose(fp);
    }
}

static void refreshGrid(fmd_t *md, int reverse)
{
    fmd_ituple_t ic, jc;
    int d;

    /* iterate over all cells */
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {

        /* iterate over all particles in cell ic */

        cell_t *cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]];
        int pind = 0;  /* particle index */

        while (pind < cell->parts_num)
        {
            particle_t *p = &cell->parts[pind];

            for (d=0; d<3; d++)
            {
                jc[d] = (int)floor(p->core.x[d] / md->cellh[d]) - md->SubDomain.ic_global_firstcell[d] + md->SubDomain.ic_start[d];

                if (jc[d] < 0)
                {
                    if (reverse && md->PBC[d] && md->ns[d] > 1 &&
                        md->SubDomain.is[d]==md->ns[d]-1 && ic[d] >= md->SubDomain.ic_stop[d]-md->SubDomain.ic_start[d])
                    {
                        p->core.x[d] += md->l[d];
                        jc[d] = ic[d] + md->SubDomain.ic_start[d];
                    }
                    else
                    {
                        fprintf(stderr,"ERROR: Unexpected particle position!\n");
                        MPI_Abort(MPI_COMM_WORLD, ERROR_UNEXPECTED_PARTICLE_POSITION);
                    }
                }
                else
                    if (jc[d] >= md->SubDomain.cell_num[d])
                    {
                        if (reverse && md->PBC[d] && md->ns[d] > 1 &&
                            md->SubDomain.is[d]==0 && ic[d] < 2*md->SubDomain.ic_start[d])
                        {
                            p->core.x[d] -= md->l[d];
                            jc[d] = ic[d] - md->SubDomain.ic_start[d];
                        }
                        else
                        {
                            fprintf(stderr, "ERROR: Unexpected particle position!\n");
                            MPI_Abort(MPI_COMM_WORLD, ERROR_UNEXPECTED_PARTICLE_POSITION);
                        }
                    }
            }

            if ((ic[0] != jc[0]) || (ic[1] != jc[1]) || (ic[2] != jc[2]))
            {
                INSERT_PARTICLE_IN_CELL(*p, md->SubDomain.grid[jc[0]][jc[1]][jc[2]]);
                REMOVE_PARTICLE_FROM_CELL(*cell, pind);
            }
            else
                pind++;
        }
    }

    /* now particles in ghost cells migrate to neighbour subdomains */
    fmd_particles_migrate(md);
}

void rescaleVelocities(fmd_t *md)
{
    fmd_ituple_t ic;
    int d;
    fmd_real_t scale;

    scale = sqrt(md->DesiredTemperature / md->GlobalTemperature);

    cell_t *cell;
    int i;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
            for (d=0; d<3; d++)
                cell->parts[i].core.v[d] *= scale;

    md->GlobalTemperature = md->DesiredTemperature;
}

/*void restoreBackups(fmd_t *md)
{
    fmd_ituple_t ic;
    int d;
    cell_t *cell;
    int pind;

    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], pind=0; pind < cell->parts_num; pind++)
        {
            particle_core_t *pc = &cell->parts[pind].core;

            for (d=0; d<3; d++)
            {
                pc->x[d] = pc->x_bak[d];
                pc->v[d] = pc->v_bak[d];
            }
        }

    refreshGrid(md, 1);

#ifdef USE_TTM
    ttm_restoreBackups();
#endif
}*/

void fmd_matt_saveConfiguration(fmd_t *md)
{
    fmd_ituple_t ic;
    XYZ_struct_t *localData, *globalData;
    int *nums, *recvcounts, *displs;
    int k;
    cell_t *cell;
    int pind;

    localData = (XYZ_struct_t *)malloc(md->SubDomain.NumberOfParticles * sizeof(XYZ_struct_t));
    k=0;
    LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], pind=0; pind < cell->parts_num; pind++)
        {
            particle_core_t *pc = &cell->parts[pind].core;

            if (md->CompLocOrdParam)
            {
                localData[k].x[0] = (float)pc->x_avgd[0];
                localData[k].x[1] = (float)pc->x_avgd[1];
                localData[k].x[2] = (float)pc->x_avgd[2];
            }
            else
            {
                localData[k].x[0] = (float)pc->x[0];
                localData[k].x[1] = (float)pc->x[1];
                localData[k].x[2] = (float)pc->x[2];
            }
            localData[k].var = pc->GroupID;
            localData[k].atomkind = pc->atomkind;
            k++;
        }

    nums = (int *)malloc(md->SubDomain.numprocs * sizeof(int));
    MPI_Allgather(&md->SubDomain.NumberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
        md->MD_comm);

    md->TotalNoOfParticles = 0;
    for (k=0; k < md->SubDomain.numprocs; k++)
        md->TotalNoOfParticles += nums[k];

    if (md->Is_MD_comm_root)
    {
        int displ = 0;

        globalData = (XYZ_struct_t *)malloc(md->TotalNoOfParticles * sizeof(XYZ_struct_t));
        recvcounts = (int *)malloc(md->SubDomain.numprocs * sizeof(int));
        displs     = (int *)malloc(md->SubDomain.numprocs * sizeof(int));
        for (k=0; k < md->SubDomain.numprocs; k++)
        {
            recvcounts[k] = nums[k] * sizeof(XYZ_struct_t);
            displs[k] = displ;
            displ += recvcounts[k];
        }
    }
    free(nums);
    MPI_Gatherv(localData, md->SubDomain.NumberOfParticles * sizeof(XYZ_struct_t), MPI_CHAR,
        globalData, recvcounts, displs, MPI_CHAR, ROOTPROCESS(md->SubDomain.numprocs),
        md->MD_comm);
    free(localData);

    if (md->Is_MD_comm_root)
    {
        char configPath[MAX_PATH_LENGTH];
        char *elementName;
        int i;

        free(recvcounts);
        free(displs);

        switch (md->SaveConfigMode)
        {
            case FMD_SCM_XYZ_PARTICLESNUM:
                if (md->TotalNoOfParticles != md->_OldNumberOfParticles)
                {
                    if (md->_OldNumberOfParticles != -1) fclose(md->ConfigFilep);
                    sprintf(configPath, "%s%d.xyz", md->SaveDirectory, md->TotalNoOfParticles);
                    md->ConfigFilep = fopen(configPath, "w");
                    handleFileOpenError(md->ConfigFilep, configPath);
                    md->_OldNumberOfParticles = md->TotalNoOfParticles;
                }
                break;

            case FMD_SCM_XYZ_SEPARATE:
                sprintf(configPath, "%s%05d.xyz", md->SaveDirectory, md->_FileIndex++);
                md->ConfigFilep = fopen(configPath, "w");
                handleFileOpenError(md->ConfigFilep, configPath);
                break;

            case FMD_SCM_CSV:
                sprintf(configPath, "%s%05d.csv", md->SaveDirectory, md->_FileIndex++);
                md->ConfigFilep = fopen(configPath, "w");
                handleFileOpenError(md->ConfigFilep, configPath);
                for (i=0; i < md->potsys.atomkinds_num; i++)
                {
                    for (k=0; k < md->TotalNoOfParticles; k++)
                        if (globalData[k].atomkind == i)
                            fprintf(md->ConfigFilep, "%.2f, %.2f, %.2f, %d, %.4f\n",
                                globalData[k].x[0], globalData[k].x[1],
                                globalData[k].x[2], i, globalData[k].var);
                }
                break;

            case FMD_SCM_VTF:
            {
                int atomID = 0;
                sprintf(configPath, "%s%05d.vtf", md->SaveDirectory, md->_FileIndex++);
                md->ConfigFilep = fopen(configPath, "w");
                handleFileOpenError(md->ConfigFilep, configPath);
                for (i=0; i < md->potsys.atomkinds_num; i++)
                {
                    int count = 0;
                    for (k=0; k < md->TotalNoOfParticles; k++)
                        if (globalData[k].atomkind == i)
                            count++;
                    fprintf(md->ConfigFilep, "atom %d:%d\tname %s\n", atomID,
                            atomID+count-1, md->potsys.atomkinds[i].name);
                    atomID += count;
                }
                atomID = 0;
                for (i=0; i < md->potsys.atomkinds_num; i++)
                {
                    for (k=0; k < md->TotalNoOfParticles; k++)
                        if (globalData[k].atomkind == i)
                        {
                            fprintf(md->ConfigFilep, "%d\tbeta %.4f\n", atomID,
                             globalData[k].var);
                            atomID++;
                        }
                }
                fprintf(md->ConfigFilep, "timestep\n");
                for (i=0; i < md->potsys.atomkinds_num; i++)
                {
                    for (k=0; k < md->TotalNoOfParticles; k++)
                        if (globalData[k].atomkind == i)
                            fprintf(md->ConfigFilep, "%.2f\t%.2f\t%.2f\n",
                                globalData[k].x[0], globalData[k].x[1],
                                globalData[k].x[2]);
                }
                break;
            }
        }

        if (md->SaveConfigMode == FMD_SCM_XYZ_SEPARATE || md->SaveConfigMode == FMD_SCM_XYZ_PARTICLESNUM)
        {
            fprintf(md->ConfigFilep, "%d\n\n", md->TotalNoOfParticles);
            for (i=0; i < md->potsys.atomkinds_num; i++)
            {
                elementName = md->potsys.atomkinds[i].name;
                for (k=0; k < md->TotalNoOfParticles; k++)
                    if (globalData[k].atomkind == i)
                        fprintf(md->ConfigFilep, "%s\t%.2f\t%.2f\t%.2f\n",
                            elementName, globalData[k].x[0], globalData[k].x[1],
                            globalData[k].x[2]);
            }
        }
        free(globalData);

        if (md->SaveConfigMode == FMD_SCM_XYZ_SEPARATE || md->SaveConfigMode == FMD_SCM_CSV ||
         md->SaveConfigMode == FMD_SCM_VTF)
            fclose(md->ConfigFilep);
    }
}

void fmd_io_saveState(fmd_t *md, fmd_string_t filename)
{
    particle_core_t *is_partcores;
    int *nums;
    int i, k;
    fmd_ituple_t ic;
    char stateFilePath[MAX_PATH_LENGTH];
    FILE *fp;
    MPI_Status status;
    cell_t *cell;
    int pind;

    if (md->Is_MD_comm_root)
    {
        nums = (int *)malloc(md->SubDomain.numprocs * sizeof(int));
        MPI_Gather(&md->SubDomain.NumberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
                   ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

        md->TotalNoOfParticles = 0;
        for (k=0; k < md->SubDomain.numprocs; k++)
            md->TotalNoOfParticles += nums[k];

        sprintf(stateFilePath, "%s%s", md->SaveDirectory, filename);
        fp = fopen(stateFilePath, "w");
        handleFileOpenError(fp, stateFilePath);
        fprintf(fp, "%.16e\n", md->MD_time);
        fprintf(fp, "%d\n", md->TotalNoOfParticles);
        fprintf(fp, "%.16e\t%.16e\t%.16e\n", md->l[0], md->l[1], md->l[2]);
        fprintf(fp, "%d %d %d\n", md->PBC[0], md->PBC[1], md->PBC[2]);

        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], pind=0; pind < cell->parts_num; pind++)
            {
                particle_core_t *pc = &cell->parts[pind].core;

                fprintf(fp, "%s %d\n", md->potsys.atomkinds[pc->atomkind].name, pc->GroupID);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", pc->x[0], pc->x[1], pc->x[2]);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", pc->v[0], pc->v[1], pc->v[2]);
            }

        for (i=0; i < ROOTPROCESS(md->SubDomain.numprocs); i++)
        {
            is_partcores = (particle_core_t *)malloc(nums[i] * sizeof(particle_core_t));
            MPI_Recv(is_partcores, nums[i] * sizeof(particle_core_t), MPI_CHAR, i, 150,
                     md->MD_comm, &status);

            for (k=0; k < nums[i]; k++)
            {
                particle_core_t *pc = is_partcores + k;

                fprintf(fp, "%s %d\n", md->potsys.atomkinds[pc->atomkind].name, pc->GroupID);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", pc->x[0], pc->x[1], pc->x[2]);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", pc->v[0], pc->v[1], pc->v[2]);
            }

            free(is_partcores);
        }

        free(nums);
        fclose(fp);
    }
    else
    {
        MPI_Gather(&md->SubDomain.NumberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
                   ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        is_partcores = (particle_core_t *)malloc(md->SubDomain.NumberOfParticles * sizeof(particle_core_t));

        k = 0;
        LOOP3D(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], pind=0; pind < cell->parts_num; pind++)
                is_partcores[k++] = cell->parts[pind].core;

        MPI_Send(is_partcores, md->SubDomain.NumberOfParticles * sizeof(particle_core_t), MPI_CHAR,
                 ROOTPROCESS(md->SubDomain.numprocs), 150, md->MD_comm);
        free(is_partcores);
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
        if (md->SubDomain.myrank == ROOTPROCESS(md->SubDomain.numprocs))
            md->Is_MD_comm_root = FMD_TRUE;
    }
}

static void create_mpi_types(fmd_t *md)
{
    MPI_Type_vector(1, 3, 0, MPI_INT, &md->mpi_types.mpi_ituple);
    MPI_Type_commit(&md->mpi_types.mpi_ituple);
}

static void free_mpi_types(fmd_t *md)
{
    MPI_Type_free(&md->mpi_types.mpi_ituple);
}

fmd_t *fmd_create()
{
    fmd_t *md = (fmd_t *)malloc(sizeof(fmd_t));

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
    md->MD_time = 0.0;
    md->time_iteration = 0;
    md->SaveDirectory[0] = '\0';
    md->CompLocOrdParam = FMD_FALSE;
    md->SubDomain.grid = NULL;
    md->TotalNoOfParticles = 0;
    md->TotalNoOfMolecules = 0;
    md->ActiveGroup = -1;             // all groups are active by default
    md->ParticlesDistributed = FMD_FALSE;
    md->GlobalGridExists = FMD_FALSE;
    md->BoxSizeDetermined = FMD_FALSE;
    md->PBCdetermined = FMD_FALSE;
    md->Is_MD_comm_root = FMD_FALSE;
    md->EventHandler = NULL;
    md->timers = NULL;
    md->timers_num = 0;
    md->turies = NULL;
    md->turies_num = 0;
    md->DesiredTemperature = 300.0;
    md->SaveConfigMode = FMD_SCM_XYZ_PARTICLESNUM;
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
    int d;

    for (d=0; d<3; d++)
    {
        md->nc[d] = (int)(md->l[d] / cutoff);
        if ((md->nc[d] < 3) && md->PBC[d])
        {
            fprintf(stderr, "ERROR: nc[%d] = %d. Under PBC, this must be greater than 2!\n", d, md->nc[d]);
            MPI_Abort(MPI_COMM_WORLD, ERROR_NC_TOO_SMALL);
        }
        md->cellh[d] = md->l[d] / md->nc[d];
    }

    if (md->Is_MD_comm_root)
    {
        md->global_grid = (cell_t ***)_fmd_array_ordinary3d_create(md->nc, sizeof(cell_t));
        assert(md->global_grid != NULL);
        /* TO-DO: handle memory error */

        _fmd_initialize_grid(md->global_grid, md->nc[0], md->nc[1], md->nc[2]);
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

void fmd_dync_setTimeStep(fmd_t *md, fmd_real_t timeStep)
{
    md->delta_t = timeStep;
}

fmd_real_t fmd_dync_getTimeStep(fmd_t *md)
{
    return md->delta_t;
}

fmd_real_t fmd_dync_getTime(fmd_t *md)
{
    return md->MD_time;
}

void fmd_dync_incTime(fmd_t *md)
{
    md->MD_time += md->delta_t;
    md->time_iteration++;
    if (md->turies_num > 0) _fmd_turies_update(md);
    if (md->EventHandler != NULL) _fmd_timer_sendTimerTickEvents(md);
}

void fmd_dync_equilibrate(fmd_t *md, int GroupID, fmd_real_t duration,
  fmd_real_t timestep, fmd_real_t strength, fmd_real_t temperature)
{
    fmd_real_t bak_mdTime, bak_DesiredTemperature;
    fmd_real_t bak_delta_t, bak_BerendsenThermostatParam;
    int bak_activeGroup;

    // make backups
    bak_mdTime = md->MD_time;
    bak_delta_t = md->delta_t;
    bak_BerendsenThermostatParam = md->BerendsenThermostatParam;
    bak_DesiredTemperature = md->DesiredTemperature;
    bak_activeGroup = md->ActiveGroup;

    // initialize
    md->MD_time = 0.0;
    fmd_dync_setTimeStep(md, timestep);
    md->DesiredTemperature = temperature;
    md->GlobalTemperature = temperature;
    md->BerendsenThermostatParam = strength;
    md->ActiveGroup = GroupID;

    // compute forces for the first time
    fmd_dync_updateForces(md);

    while (md->MD_time < duration)
    {
        // take first step of velocity Verlet integrator
        fmd_dync_VelocityVerlet_startStep(md, FMD_TRUE);

        // compute forces
        fmd_dync_updateForces(md);

        // take last step of velocity Verlet integrator
        fmd_dync_VelocityVerlet_finishStep(md);

        fmd_dync_incTime(md);
    }
    // end of the time loop

    // restore backups
    md->MD_time = bak_mdTime;
    fmd_dync_setTimeStep(md, bak_delta_t);
    md->DesiredTemperature = bak_DesiredTemperature;
    md->BerendsenThermostatParam = bak_BerendsenThermostatParam;
    md->ActiveGroup = bak_activeGroup;
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
    int *start, *stop;
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
    int d;

    cell_t *cell;
    int pind;

    LOOP3D(ic, start, stop)
        for (cell = &grid[ic[0]][ic[1]][ic[2]], pind = 0; pind < cell->parts_num; pind++)
        {
            particle_core_t *pc = &cell->parts[pind].core;

            if (GroupID == -1 || GroupID == pc->GroupID)
            {
                fmd_real_t mass = md->potsys.atomkinds[pc->atomkind].mass;
                fmd_real_t stdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

                for (d=0; d<3; d++)
                    pc->v[d] = gsl_ran_gaussian_ziggurat(rng, stdDevVelocity);
            }
        }

    gsl_rng_free(rng);
}

fmd_real_t fmd_matt_getGlobalTemperature(fmd_t *md)
{
    return md->GlobalTemperature;
}

void fmd_dync_setBerendsenThermostatParameter(fmd_t *md, fmd_real_t parameter)
{
    md->BerendsenThermostatParam = parameter;
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
