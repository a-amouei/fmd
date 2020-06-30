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

const int fmd_ThreeZeros[3] = {0, 0, 0};

static void refreshGrid(fmd_t *md, int reverse);

void cleanGridSegment(cell_t ***grid, int ic_from[3], int ic_to[3])
{
    int ic[3];
    ParticleListItem_t *item_p, **item_pp;

    ITERATE(ic, ic_from, ic_to)
    {
        item_pp = &grid[ic[0]][ic[1]][ic[2]];
        item_p = *item_pp;
        while (item_p != NULL)
        {
            removeFromList(item_pp);
            free(item_p);
            item_p = *item_pp;
        }
    }
}

void compLocOrdParam(fmd_t *md)
{
/*
    float latticeParameter = md->EAM.elements[0].latticeParameter;
    float rCutSqd = SQR(1.32 * latticeParameter);
    int ic[3], jc[3], kc[3];
    int d;
    ParticleListItem_t *item1_p, *item2_p;
    float arg, real, img;
    int Z;
    float r2, rv[3];
    float q[6][3] = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{0,1,1},{1,0,1}};
    int i;

    if (md->LOP_iteration == 0)
        ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
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

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item1_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            real = img = 0.;
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
                                r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);
                                if (r2 < rCutSqd)
                                {
                                    Z++;
                                    for (i=0; i<6; i++)
                                    {
                                        arg = 0.;
                                        for (d=0; d<3; d++)
                                            arg += q[i][d] * rv[d];
                                        real += cosf(arg);
                                        img  += sinf(arg);
                                    }
                                }
                            }
                    }
                }
            }
            item1_p->P.LocOrdParam += (SQR(real)+SQR(img)) / SQR(Z);
            for (d=0; d<3; d++)
                item1_p->P.x_avgd[d] += item1_p->P.x[d];
        }

    if (md->LOP_iteration != md->LOP_period) return;

    fmd_ghostparticles_update_LocOrdParam(md);

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item1_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            real = item1_p->P.LocOrdParam;
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
                                r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);
                                if (r2 < rCutSqd)
                                {
                                    Z++;
                                    real += item2_p->P.LocOrdParam;

                                }
                            }
                    }
                }
            }
            item1_p->LocOrdParamAvg = real / ((Z+1) * 36 * md->LOP_period);
            for (d=0; d<3; d++)
                item1_p->P.x_avgd[d] /= md->LOP_period;
        }

    md->LOP_iteration = 0;
*/
}

void fmd_dync_VelocityVerlet_startStep(fmd_t *md, fmd_bool_t UseThermostat)
{
    int ic[3];
    int d;
    ParticleListItem_t *item_p, **item_pp;
    double velocityScale, mass;
    double x;
    int itemDestroyed;

    if (UseThermostat) velocityScale = sqrt(1 + md->delta_t / md->BerendsenThermostatParam *
                       (md->DesiredTemperature / md->GlobalTemperature - 1));
    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {
        // iterate over all items in cell ic
        item_pp = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]];
        item_p = *item_pp;
        while (item_p != NULL)
        {
            if (!(md->ActiveGroup == -1 || item_p->P.GroupID == md->ActiveGroup))
            {
                item_pp = &item_p->next_p;
                item_p = *item_pp;
                continue;
            }

            itemDestroyed = 0;
            mass = md->potsys.atomkinds[item_p->P.atomkind].mass;

            for (d=0; d<3; d++)
            {
                if (md->UseAutoStep)
                {
                    item_p->P.v_bak[d] = item_p->P.v[d];
                    item_p->P.x_bak[d] = item_p->P.x[d];
                }
                if (UseThermostat) item_p->P.v[d] *= velocityScale;
                item_p->P.v[d] += md->delta_t * 0.5 / mass * item_p->F[d];
                x = item_p->P.x[d] + md->delta_t * item_p->P.v[d];

                if ( (md->ns[d] == 1) && ((x < 0.0) || (x >= md->l[d])) )
                {
                    if (!md->PBC[d])
                    {
                        removeFromList(item_pp);
                        free(item_p);
                        (md->SubDomain.NumberOfParticles)--;
                        itemDestroyed = 1;
                        break;
                    }
                    else
                        if (x < 0.0) x += md->l[d]; else x -= md->l[d];
                }
                item_p->P.x[d] = x;
            }

            if (!itemDestroyed)
                item_pp = &item_p->next_p;
            item_p = *item_pp;
        }
    }
    refreshGrid(md, 0);
}

int fmd_dync_VelocityVerlet_finishStep(fmd_t *md)
{
    int ic[3];
    int d;
    ParticleListItem_t *item_p;
    double m_vSqd_Sum = 0, m_vSqd_SumSum;
    double mass;
    int returnVal = 0;
    int particlesNum = 0;
    double momentumSum[3] = {0., 0., 0.};

    for (ic[0] = md->SubDomain.ic_start[0]; ic[0] < md->SubDomain.ic_stop[0]; ic[0]++)
    {
        for (ic[1] = md->SubDomain.ic_start[1]; ic[1] < md->SubDomain.ic_stop[1]; ic[1]++)
            for (ic[2] = md->SubDomain.ic_start[2]; ic[2] < md->SubDomain.ic_stop[2]; ic[2]++)
                for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                {
                    if (!(md->ActiveGroup == -1 || item_p->P.GroupID == md->ActiveGroup))
                        continue;
                    particlesNum++;
                    mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
                    for (d=0; d<3; d++)
                    {
                        item_p->P.v[d] += md->delta_t * 0.5 / mass * item_p->F[d];
                        momentumSum[d] += mass * item_p->P.v[d];
                    }
                    m_vSqd_Sum += mass * ( SQR(item_p->P.v[0]) +
                                           SQR(item_p->P.v[1]) +
                                           SQR(item_p->P.v[2]) );
                }
    }
    MPI_Reduce(momentumSum, md->TotalMomentum, 3, MPI_DOUBLE, MPI_SUM,
               ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
    MPI_Allreduce(&particlesNum, &(md->ActiveGroupParticlesNum), 1, MPI_INT, MPI_SUM, md->MD_comm);

    if (md->ActiveGroup == -1)
        md->TotalNoOfParticles = md->ActiveGroupParticlesNum;

    MPI_Allreduce(&m_vSqd_Sum, &m_vSqd_SumSum, 1, MPI_DOUBLE, MPI_SUM, md->MD_comm);
    md->TotalKineticEnergy = 0.5 * m_vSqd_SumSum;
    md->TotalMDEnergy = md->TotalKineticEnergy + md->TotalPotentialEnergy;

    if (md->UseAutoStep)
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
    }
    md->GlobalTemperature = m_vSqd_SumSum / (3.0 * md->ActiveGroupParticlesNum * K_BOLTZMANN);

    return returnVal;
}

// not correct under periodic boundary conditions
// see [J. Chem. Phys. 131, 154107 (2009)]
static double compVirial_internal(fmd_t *md)
{
    int ic[3];
    ParticleListItem_t *item_p;
    double virial = 0.0;
    double virial_global;

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            virial += item_p->P.x[0] * item_p->F[0] +
                      item_p->P.x[1] * item_p->F[1] +
                      item_p->P.x[2] * item_p->F[2];
        }
    MPI_Reduce(&virial, &virial_global, 1, MPI_DOUBLE, MPI_SUM,
      ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

    return virial_global;
}

void createCommunicators(fmd_t *md)
{
    int mdnum, i;
    MPI_Group world_group, MD_group;
    int *ranks;

    // create MD_comm
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

cell_t ***fmd_internal_createGrid(int cell_num[3])
{
    cell_t ***grid;
    int i, j, k;

    grid = (cell_t ***)malloc(cell_num[0]*sizeof(cell_t **));
    for (i=0; i < cell_num[0]; i++)
    {
        grid[i] = (cell_t **)malloc(cell_num[1]*sizeof(cell_t *));
        for (j=0; j < cell_num[1]; j++)
        {
            grid[i][j] = (cell_t *)malloc(cell_num[2]*sizeof(cell_t));
            for (k=0; k < cell_num[2]; k++)
                grid[i][j][k] = NULL;
        }
    }
    return grid;
}

void findLimits(fmd_t *md, double LowerLimit[3], double UpperLimit[3])
{
    ParticleListItem_t *item_p;
    int ic[3];
    int d;
    double LocalLower[3], LocalUpper[3];

    LocalLower[0] = LocalLower[1] = LocalLower[2] = DBL_MAX;
    LocalUpper[0] = LocalUpper[1] = LocalUpper[2] = DBL_MIN;

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            for (d=0; d<3; d++)
            {
                if (item_p->P.x[d] < LocalLower[d])
                    LocalLower[d] = item_p->P.x[d];
                if (item_p->P.x[d] > LocalUpper[d])
                    LocalUpper[d] = item_p->P.x[d];
            }
    MPI_Allreduce(LocalLower, LowerLimit, 3, MPI_DOUBLE, MPI_MIN, md->MD_comm);
    MPI_Allreduce(LocalUpper, UpperLimit, 3, MPI_DOUBLE, MPI_MAX, md->MD_comm);
}

void fmd_internal_freeGrid(cell_t ***grid, int *cell_num)
{
    int i, j;

    cleanGridSegment(grid, fmd_ThreeZeros, cell_num);
    for (i=0; i < cell_num[0]; i++)
    {
        for (j=0; j < cell_num[1]; j++)
            free(grid[i][j]);
        free(grid[i]);
    }
    free(grid);
}

int getListLength(ParticleListItem_t *root_p)
{
    int i = 0;

    for ( ; root_p != NULL; root_p = root_p->next_p)
        i++;
    return i;
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

void fmd_matt_addVelocity(fmd_t *md, int GroupID, double vx, double vy, double vz)
{
    cell_t ***grid;
    int *start, *stop;
    int ic[3];
    ParticleListItem_t *item_p;

    if (md->ParticlesDistributed)
    {
        grid = md->SubDomain.grid;
        start = md->SubDomain.ic_start;
        stop = md->SubDomain.ic_stop;
    }
    else
    {
        start = fmd_ThreeZeros;
        if (md->Is_MD_comm_root)
        {
            grid = md->global_grid;
            stop = md->nc;
        }
        else
            stop = start;
    }

    ITERATE(ic, start, stop)
        for (item_p = grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            if (GroupID == -1 || GroupID == item_p->P.GroupID)
            {
                item_p->P.v[0] += vx;
                item_p->P.v[1] += vy;
                item_p->P.v[2] += vz;
            }
}

void fmd_matt_distribute(fmd_t *md)
{
    ParticleListItem_t *item_p;
    int i, k, d, nct, sum_length;
    int ic[3], *ic_length;
    particle_t *is_particles;

    if (md->SubDomain.grid == NULL) fmd_subd_init(md);

    if (md->Is_MD_comm_root)
    {
        int r, w;
        int is[3], global_icstart[3], global_icstop[3];
        ParticleListItem_t **item_pp;
        double m_vSqd_Sum = 0.0;
        double mass;

#ifdef USE_TTM
        ttm_comp_min_atomsNo(global_grid, s_p);
#endif

        ITERATE(ic, fmd_ThreeZeros, md->nc)
            for (item_p=md->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
                for (d=0; d<3; d++)
                    m_vSqd_Sum += mass * SQR(item_p->P.v[d]);
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
            ITERATE(ic, global_icstart, global_icstop)
            {
                ic_length[k] = getListLength(md->global_grid[ic[0]][ic[1]][ic[2]]);
                sum_length += ic_length[k++];
            }
            ic_length[nct] = sum_length;
            MPI_Send(ic_length, nct+1, MPI_INT, i, 50, md->MD_comm);
            free(ic_length);
            sum_length *= sizeof(particle_t);
            is_particles = (particle_t *)malloc(sum_length);
            k = 0;
            ITERATE(ic, global_icstart, global_icstop)
                for (item_p=md->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    is_particles[k++] = item_p->P;
            MPI_Send(is_particles, sum_length, MPI_CHAR, i, 51, md->MD_comm);
            free(is_particles);
        }

        md->SubDomain.NumberOfParticles = 0;
        ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        {
            item_pp = &md->global_grid[ ic[0] - md->SubDomain.ic_start[0] + md->SubDomain.ic_global_firstcell[0] ]
                                      [ ic[1] - md->SubDomain.ic_start[1] + md->SubDomain.ic_global_firstcell[1] ]
                                      [ ic[2] - md->SubDomain.ic_start[2] + md->SubDomain.ic_global_firstcell[2] ];
            item_p = *item_pp;

            while (item_p != NULL)
            {
                removeFromList(item_pp);
                item_p->neighbors = NULL;
                fmd_insertInList(&md->SubDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
                ++(md->SubDomain.NumberOfParticles);
                item_p = *item_pp;
            }


        }
        fmd_internal_freeGrid(md->global_grid, md->nc);
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
        sum_length *= sizeof(particle_t);
        is_particles = (particle_t *)malloc(sum_length);
        MPI_Recv(is_particles, sum_length, MPI_CHAR,
                 ROOTPROCESS(md->SubDomain.numprocs), 51, md->MD_comm, &status);
        kreceive = k = 0;
        ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        {
            for (i=0; i<ic_length[kreceive]; i++)
            {
                item_p = (ParticleListItem_t *)malloc(sizeof(ParticleListItem_t));
                item_p->P = is_particles[k++];
                item_p->neighbors = NULL;
                fmd_insertInList(&md->SubDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
            }
            kreceive++;
        }
        free(ic_length);
        free(is_particles);
    }

    MPI_Bcast(&md->TotalNoOfParticles, 1, MPI_UNSIGNED, ROOTPROCESS(md->SubDomain.numprocs),
              md->MD_comm);
    MPI_Bcast(&md->TotalNoOfMolecules, 1, MPI_UNSIGNED, ROOTPROCESS(md->SubDomain.numprocs),
              md->MD_comm);
    MPI_Bcast(&md->GlobalTemperature, 1, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs),
              md->MD_comm);

    if (md->TotalNoOfMolecules > 0) fmd_matt_updateNeighbors(md);

    md->TotalKineticEnergy = 3.0/2.0 * md->TotalNoOfParticles * K_BOLTZMANN * md->GlobalTemperature;
    md->GlobalGridExists = FMD_FALSE;
    md->ParticlesDistributed = FMD_TRUE;
}

void fmd_insertInList(ParticleListItem_t **root_pp, ParticleListItem_t *item_p)
{
    item_p->next_p = *root_pp;
    *root_pp = item_p;
}

void fmd_io_loadState(fmd_t *md, fmd_string_t file, fmd_bool_t useTime)
{
    ParticleListItem_t *item_p;
    FILE *fp;
    char name[3];
    int i, j, d;
    int ic[3];
    double stateFileTime;
    int particlesNum;
    double l0, l1, l2;
    int PBC0, PBC1, PBC2;

    if (md->Is_MD_comm_root)
    {
        fp = fopen(file, "r");
        handleFileOpenError(fp, file);
        fscanf(fp, "%lf", &stateFileTime);
        if (useTime)
            md->MD_time = stateFileTime;
        fscanf(fp, "%d\n", &particlesNum);
        md->TotalNoOfParticles += particlesNum;
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
        MPI_Bcast(&md->l, 3, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
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

    if (useTime)
        MPI_Bcast(&md->MD_time, 1, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        for (i=0; i < particlesNum; i++)
        {
            item_p = (ParticleListItem_t *)malloc(sizeof(ParticleListItem_t));
            fscanf(fp, "%s%d", name, &item_p->P.GroupID);
            for (j=0; j < md->potsys.atomkinds_num; j++)
                if (strcmp(name, md->potsys.atomkinds[j].name) == 0)
                {
                    item_p->P.atomkind = j;
                    break;
                }
            // TO-DO: what if the name doesn't exist in potsys?
            fscanf(fp, "%lf%lf%lf", &item_p->P.x[0], &item_p->P.x[1], &item_p->P.x[2]);
            fscanf(fp, "%lf%lf%lf", &item_p->P.v[0], &item_p->P.v[1], &item_p->P.v[2]);
            for (d=0; d<3; d++)
                ic[d] = (int)floor(item_p->P.x[d] / md->cellh[d]);
            fmd_insertInList(&(md->global_grid[ic[0]][ic[1]][ic[2]]), item_p);
        }
        fclose(fp);
    }
}

static void refreshGrid(fmd_t *md, int reverse)
{
    int ic[3], jc[3];
    int d;
    ParticleListItem_t **item_pp;
    ParticleListItem_t *item_p;

    // iterate over all cells(lists)
    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
    {
        // iterate over all items in cell ic
        item_pp = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]];
        item_p = *item_pp;
        while (item_p != NULL)
        {
            for (d=0; d<3; d++)
            {
                jc[d] = (int)floor(item_p->P.x[d] / md->cellh[d]) - md->SubDomain.ic_global_firstcell[d] + md->SubDomain.ic_start[d];
                if (jc[d] < 0)
                {
                    if (reverse && md->PBC[d] && md->ns[d] > 1 &&
                        md->SubDomain.is[d]==md->ns[d]-1 && ic[d] >= md->SubDomain.ic_stop[d]-md->SubDomain.ic_start[d])
                    {
                        item_p->P.x[d] += md->l[d];
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
                            item_p->P.x[d] -= md->l[d];
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
                removeFromList(item_pp);
                fmd_insertInList(&md->SubDomain.grid[jc[0]][jc[1]][jc[2]], item_p);
            }
            else
                item_pp = &item_p->next_p;

            item_p = *item_pp;
        }
    }

    // now particles in ghost cells migrate to neighbour subdomains
    fmd_particles_migrate(md);
}

void removeFromList(ParticleListItem_t **item_pp)
{
    *item_pp = (*item_pp)->next_p;
}

void rescaleVelocities(fmd_t *md)
{
    int ic[3];
    int d;
    ParticleListItem_t *item_p;
    double scale;

    scale = sqrt(md->DesiredTemperature / md->GlobalTemperature);

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            for (d=0; d<3; d++)
                item_p->P.v[d] *= scale;
    md->GlobalTemperature = md->DesiredTemperature;
}

void restoreBackups(fmd_t *md)
{
    int ic[3], d;
    ParticleListItem_t *item_p;

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            for (d=0; d<3; d++)
            {
                item_p->P.x[d] = item_p->P.x_bak[d];
                item_p->P.v[d] = item_p->P.v_bak[d];
            }
        }
    refreshGrid(md, 1);

#ifdef USE_TTM
    ttm_restoreBackups();
#endif
}

void fmd_matt_saveConfiguration(fmd_t *md)
{
    int ic[3];
    ParticleListItem_t *item_p;
    XYZ_struct_t *localData, *globalData;
    int *nums, *recvcounts, *displs;
    int k;

    localData = (XYZ_struct_t *)malloc(md->SubDomain.NumberOfParticles * sizeof(XYZ_struct_t));
    k=0;
    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            if (md->CompLocOrdParam)
            {
                localData[k].x[0] = (float)item_p->P.x_avgd[0];
                localData[k].x[1] = (float)item_p->P.x_avgd[1];
                localData[k].x[2] = (float)item_p->P.x_avgd[2];
            }
            else
            {
                localData[k].x[0] = (float)item_p->P.x[0];
                localData[k].x[1] = (float)item_p->P.x[1];
                localData[k].x[2] = (float)item_p->P.x[2];
            }
            localData[k].var = item_p->P.GroupID;
            localData[k].atomkind = item_p->P.atomkind;
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
    particle_t *is_particles, *P_p;
    int *nums;
    int i, k, ic[3];
    ParticleListItem_t *item_p;
    char stateFilePath[MAX_PATH_LENGTH];
    FILE *fp;
    MPI_Status status;

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
        ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                fprintf(fp, "%s %d\n", md->potsys.atomkinds[item_p->P.atomkind].name, item_p->P.GroupID);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", item_p->P.x[0], item_p->P.x[1], item_p->P.x[2]);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", item_p->P.v[0], item_p->P.v[1], item_p->P.v[2]);
            }
        for (i=0; i < ROOTPROCESS(md->SubDomain.numprocs); i++)
        {
            is_particles = (particle_t *)malloc(nums[i] * sizeof(particle_t));
            MPI_Recv(is_particles, nums[i] * sizeof(particle_t), MPI_CHAR, i, 150,
                md->MD_comm, &status);
            for (k=0; k < nums[i]; k++)
            {
                P_p = is_particles + k;
                fprintf(fp, "%s %d\n", md->potsys.atomkinds[P_p->atomkind].name, P_p->GroupID);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", P_p->x[0], P_p->x[1], P_p->x[2]);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", P_p->v[0], P_p->v[1], P_p->v[2]);
            }
            free(is_particles);
        }
        free(nums);
        fclose(fp);
    }
    else
    {
        MPI_Gather(&md->SubDomain.NumberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
            ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        is_particles = (particle_t *)malloc(md->SubDomain.NumberOfParticles * sizeof(particle_t));
        k = 0;
        ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
            for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                is_particles[k++] = item_p->P;
        MPI_Send(is_particles, md->SubDomain.NumberOfParticles * sizeof(particle_t), MPI_CHAR,
            ROOTPROCESS(md->SubDomain.numprocs), 150, md->MD_comm);
        free(is_particles);
    }
}

void fmd_matt_setDesiredTemperature(fmd_t *md, double DesiredTemperature)
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

    // this must be the last statement before return
    md->WallTimeOrigin = MPI_Wtime();

    return md;
}

void fmd_box_setSize(fmd_t *md, double sx, double sy, double sz)
{
    if (!md->GlobalGridExists)
    {
        md->l[0] = sx;
        md->l[1] = sy;
        md->l[2] = sz;
        md->BoxSizeDetermined = FMD_TRUE;
    }
}

double fmd_proc_getWallTime(fmd_t *md)
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

void fmd_box_createGrid(fmd_t *md, double cutoff)
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
        md->global_grid = fmd_internal_createGrid(md->nc);
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

void fmd_dync_setTimeStep(fmd_t *md, double timeStep)
{
    md->delta_t = timeStep;
}

double fmd_dync_getTimeStep(fmd_t *md)
{
    return md->delta_t;
}

double fmd_dync_getTime(fmd_t *md)
{
    return md->MD_time;
}

void fmd_dync_incTime(fmd_t *md)
{
    md->MD_time += md->delta_t;
    if (md->EventHandler != NULL) fmd_timer_sendTimerTickEvents(md);
}

void fmd_dync_equilibrate(fmd_t *md, int GroupID, double duration,
  double timestep, double strength, double temperature)
{
    double bak_mdTime, bak_DesiredTemperature;
    double bak_delta_t, bak_BerendsenThermostatParam;
    int bak_activeGroup;

    // make backups
    bak_mdTime = md->MD_time;
    bak_delta_t = md->delta_t;
    bak_BerendsenThermostatParam = md->BerendsenThermostatParam;
    bak_DesiredTemperature = md->DesiredTemperature;
    bak_activeGroup = md->ActiveGroup;

    // initialize
    md->MD_time = 0.0;
    md->delta_t = timestep;
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

        md->MD_time += md->delta_t;
        if (md->EventHandler != NULL) fmd_timer_sendTimerTickEvents(md);
    }
    // end of the time loop

    // restore backups
    md->MD_time = bak_mdTime;
    md->delta_t = bak_delta_t;
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

double fmd_matt_getTotalEnergy(fmd_t *md)
{
    return md->TotalKineticEnergy + md->TotalPotentialEnergy;
}

void fmd_matt_giveTemperature(fmd_t *md, int GroupID)
{
    cell_t ***grid;
    int *start, *stop;
    int ic[3];
    ParticleListItem_t *item_p;

    if (md->ParticlesDistributed)
    {
        grid = md->SubDomain.grid;
        start = md->SubDomain.ic_start;
        stop = md->SubDomain.ic_stop;
    }
    else
    {
        start = fmd_ThreeZeros;
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

    ITERATE(ic, start, stop)
        for (item_p = grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            if (GroupID == -1 || GroupID == item_p->P.GroupID)
            {
                double mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
                double stdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

                for (d=0; d<3; d++)
                    item_p->P.v[d] = gsl_ran_gaussian_ziggurat(rng, stdDevVelocity);
            }
        }

    gsl_rng_free(rng);
}

double fmd_matt_getGlobalTemperature(fmd_t *md)
{
    return md->GlobalTemperature;
}

void fmd_dync_setBerendsenThermostatParameter(fmd_t *md, double parameter)
{
    md->BerendsenThermostatParam = parameter;
}

void fmd_free(fmd_t *md)
{
    fmd_subd_free(md);
    fmd_potsys_free(md);
    fmd_timer_free(md);

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
