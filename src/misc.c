/*
  misc.c: This file is part of Free Molecular Dynamics

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

#include <string.h>
#include <tgmath.h>
#include <omp.h>
#include "matter.h"
#include "misc.h"
#include "fmd-private.h"
#include "md-ghost.h"
#include "timer.h"
#include "turi.h"
#include "general.h"

typedef struct
{
    fmd_rtuple_t x;
    fmd_rtuple_t v;
    int GroupID;
    unsigned atomkind;
} state_atom_t;

static char formatstr_3xpoint16e[] = "%.16e\t%.16e\t%.16e\n";

void _fmd_createGlobalGrid(fmd_t *md)
{
    fmd_real_t cutoff = _fmd_pot_get_largest_cutoff(&md->potsys);

    for (int d=0; d < DIM; d++)
    {
        md->nc[d] = (int)(md->l[d] / cutoff);
        md->cellh[d] = md->l[d] / md->nc[d];
    }

    _fmd_cellinfo_init(&md->cellinfo);

    if (md->Is_MD_comm_root)
    {
        unsigned nc = md->nc[0] * md->nc[1] * md->nc[2];

        md->ggrid = m_alloc(nc * sizeof(cell_t));

        for (int ic=0; ic < nc; ic++)
            _fmd_cell_init(&md->cellinfo, md->ggrid + ic);
    }
}

static void compLocOrdParam(fmd_t *md)
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
        LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
            for (item1_p = ARRAY_ELEMENT(md->Subdomain.grid, ic); item1_p != NULL; item1_p = item1_p->next_p)
            {
                item1_p->P.LocOrdParam = 0.;
                for (d=0; d<3; d++)
                    item1_p->P.x_avgd[d] = 0.0;
            }

    (md->LOP_iteration)++;

    for (i=0; i<6; i++)
        for (d=0; d<3; d++)
            q[i][d] *= 4.0 * M_PI / latticeParameter * q[i][d];

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (item1_p = ARRAY_ELEMENT(md->Subdomain.grid, ic); item1_p != NULL; item1_p = item1_p->next_p)
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
                        for (item2_p = ARRAY_ELEMENT(md->Subdomain.grid, jc); item2_p != NULL; item2_p = item2_p->next_p)
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

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (item1_p = ARRAY_ELEMENT(md->Subdomain.grid, ic); item1_p != NULL; item1_p = item1_p->next_p)
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
                        for (item2_p = ARRAY_ELEMENT(md->Subdomain.grid, jc); item2_p != NULL; item2_p = item2_p->next_p)
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
    fmd_real_t virial = 0.0;
    fmd_real_t virial_global;
    unsigned i;
    cell_t *c;

    for (int ic=0; ic < md->subd.nc; ic++)
        for (c = md->subd.grid + ic, i=0; i < c->parts_num; i++)
            virial += POS(c, i, 0) * FRC(c, i, 0) +
                      POS(c, i, 1) * FRC(c, i, 1) +
                      POS(c, i, 2) * FRC(c, i, 2);

    MPI_Reduce(&virial, &virial_global, 1, FMD_MPI_REAL, MPI_SUM, RANK0, md->MD_comm);

    return virial_global;
}

static void createCommunicators(fmd_t *md)
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

static void identifyProcess(fmd_t *md)
{
    int mdnum = md->ns[0] * md->ns[1] * md->ns[2];

    if (md->world_rank < mdnum)
        md->Is_MD_process = true;
    else
        md->Is_MD_process = false;
}

void fmd_io_loadState(fmd_t *md, fmd_string_t path, bool UseTime)
{
    FILE *fp;
    char name[3];
    fmd_real_t StateFileTime;
    unsigned ParticlesNum;
    fmd_real_t l0, l1, l2;
    int PBC0, PBC1, PBC2;

    if (md->Is_MD_comm_root)
    {
        fp = f_open(path, "r");

        int numread = fscanf(fp, "%lf", &StateFileTime);
        assert(numread == 1); /* TO-DO: handle error */

        if (UseTime) md->time = StateFileTime;

        numread = fscanf(fp, "%u\n", &ParticlesNum);
        assert(numread == 1); /* TO-DO: handle error */
        numread = fscanf(fp, "%lf%lf%lf", &l0, &l1, &l2);
        assert(numread == 3); /* TO-DO: handle error */
        numread = fscanf(fp, "%d%d%d", &PBC0, &PBC1, &PBC2);
        assert(numread == 3); /* TO-DO: handle error */
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
        md->BoxSizeDetermined = true;
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
        md->PBCdetermined = true;
    }

    if (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (UseTime)
        MPI_Bcast(&md->time, 1, FMD_MPI_REAL, RANK0, md->MD_comm);

    if (md->Is_MD_comm_root)
    {
        for (unsigned i=0; i < ParticlesNum; i++)
        {
            int GroupID;
            unsigned atomkind;

            int numread = fscanf(fp, "%s%d", name, &GroupID);
            assert(numread == 2); /* TO-DO: handle error */
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

            numread = fscanf(fp, "%lf%lf%lf", &x[0], &x[1], &x[2]);
            assert(numread == 3); /* TO-DO: handle error */
            numread = fscanf(fp, "%lf%lf%lf", &v[0], &v[1], &v[2]);
            assert(numread == 3); /* TO-DO: handle error */

            fmd_ituple_t ic;

            for (int d=0; d<DIM; d++)
                ic[d] = (int)floor(x[d] / md->cellh[d]);

            cell_t *c = md->ggrid + INDEX_FLAT(ic, md->nc);

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

    md->KineticEnergyUpdated = false;
}

void _fmd_refreshGrid(fmd_t *md)
{
    /* iterate over all cells */
    for (int ic=0; ic < md->subd.nc; ic++)
    {
        /* iterate over all particles in cell ic */

        cell_t *c = md->subd.grid + ic;
        int i = 0;  /* particle index */
        fmd_ituple_t jcv;

        while (i < c->parts_num)
        {
            for (int d=0; d<DIM; d++)
            {
                jcv[d] = (int)floor(POS(c, i, d) / md->cellh[d]) - md->subd.ic_global_firstcell[d]
                         + md->subd.ic_start[d];

                if (jcv[d] < 0 || jcv[d] >= md->subd.cell_num[d])
                {
                    fprintf(stderr, "ERROR: Unexpected particle position!\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR_UNEXPECTED_PARTICLE_POSITION);
                }
            }

            cell_t *c2 = ARRAY_ELEMENT(md->subd.gridp, jcv);

            if (c != c2)
            {
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

void fmd_io_saveState(fmd_t *md, fmd_string_t filename)
{
    unsigned *nums;
    char StateFilePath[MAX_PATH_LENGTH];
    FILE *fp;
    MPI_Status status;
    cell_t *c;
    unsigned pi;  /* particle index, not the famous pi number! */

    if (!md->ParticlesDistributed)
    {
        if (md->Is_MD_comm_root)
        {
            sprintf(StateFilePath, "%s%s", md->SaveDirectory, filename);
            fp = f_open(StateFilePath, "w");

            fprintf(fp, "%.16e\n", md->time);
            fprintf(fp, "%u\n", md->TotalNoOfParticles);
            fprintf(fp, formatstr_3xpoint16e, md->l[0], md->l[1], md->l[2]);
            fprintf(fp, "%d %d %d\n", md->PBC[0], md->PBC[1], md->PBC[2]);

            int nc = md->nc[0] * md->nc[1] * md->nc[2];

            for (int ic=0; ic < nc; ic++)
                for (c = md->ggrid + ic, pi=0; pi < c->parts_num; pi++)
                {
                    fprintf(fp, "%s %d\n", md->potsys.atomkinds[c->atomkind[pi]].name, c->GroupID[pi]);
                    fprintf(fp, formatstr_3xpoint16e, POS(c, pi, 0), POS(c, pi, 1), POS(c, pi, 2));
                    fprintf(fp, formatstr_3xpoint16e, VEL(c, pi, 0), VEL(c, pi, 1), VEL(c, pi, 2));
                }

            fclose(fp);
        }

        return;
    }

    if (md->Is_MD_comm_root)
    {
        nums = (unsigned *)m_alloc(md->subd.numprocs * sizeof(unsigned));

        MPI_Gather(&md->subd.NumberOfParticles, 1, MPI_UNSIGNED, nums, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

        md->TotalNoOfParticles = 0;

        for (int k=0; k < md->subd.numprocs; k++)
            md->TotalNoOfParticles += nums[k];

        sprintf(StateFilePath, "%s%s", md->SaveDirectory, filename);
        fp = f_open(StateFilePath, "w");

        fprintf(fp, "%.16e\n", md->time);
        fprintf(fp, "%u\n", md->TotalNoOfParticles);
        fprintf(fp, formatstr_3xpoint16e, md->l[0], md->l[1], md->l[2]);
        fprintf(fp, "%d %d %d\n", md->PBC[0], md->PBC[1], md->PBC[2]);

        for (int ic=0; ic < md->subd.nc; ic++)
            for (c = md->subd.grid + ic, pi=0; pi < c->parts_num; pi++)
            {
                fprintf(fp, "%s %d\n", md->potsys.atomkinds[c->atomkind[pi]].name, c->GroupID[pi]);
                fprintf(fp, formatstr_3xpoint16e, POS(c, pi, 0), POS(c, pi, 1), POS(c, pi, 2));
                fprintf(fp, formatstr_3xpoint16e, VEL(c, pi, 0), VEL(c, pi, 1), VEL(c, pi, 2));
            }

        for (int i=1; i < md->subd.numprocs; i++)
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
        MPI_Gather(&md->subd.NumberOfParticles, 1, MPI_UNSIGNED, nums, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

        state_atom_t *states = (state_atom_t *)m_alloc(md->subd.NumberOfParticles * sizeof(state_atom_t));

        unsigned k = 0;

        for (int ic=0; ic < md->subd.nc; ic++)
            for (c = md->subd.grid + ic, pi=0; pi < c->parts_num; pi++)
            {
                for (int d=0; d<DIM; d++)
                    states[k].x[d] = POS(c, pi, d);

                for (int d=0; d<DIM; d++)
                    states[k].v[d] = VEL(c, pi, d);

                states[k].GroupID = c->GroupID[pi];
                states[k].atomkind = c->atomkind[pi];

                k++;
            }

        MPI_Send(states, md->subd.NumberOfParticles, md->mpi_types.mpi_statea, RANK0, 150, md->MD_comm);

        free(states);
    }
}

void fmd_box_setPBC(fmd_t *md, bool PBCx, bool PBCy, bool PBCz)
{
    md->PBC[0] = PBCx;
    md->PBC[1] = PBCy;
    md->PBC[2] = PBCz;
    md->PBCdetermined = true;
}

bool fmd_box_setSubdomains(fmd_t *md, int dimx, int dimy, int dimz)
{
    md->ns[0] = dimx;
    md->ns[1] = dimy;
    md->ns[2] = dimz;

    identifyProcess(md);
    createCommunicators(md);

    if (md->Is_MD_process)
    {
        MPI_Comm_size(md->MD_comm, &md->subd.numprocs);
        MPI_Comm_rank(md->MD_comm, &md->subd.myrank);
        if (md->subd.myrank == RANK0)
            md->Is_MD_comm_root = true;
    }

    return md->Is_MD_process;
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

    md->MPI_initialized_by_me = false;

    int isMPIInitialized;

    MPI_Initialized(&isMPIInitialized);

    if (!isMPIInitialized)
    {
        MPI_Init(NULL, NULL);
        md->MPI_initialized_by_me = true;
    }

    md->numthreads = 1;

    MPI_Comm_size(MPI_COMM_WORLD, &(md->world_numprocs));
    MPI_Comm_rank(MPI_COMM_WORLD, &(md->world_rank));
    md->LOP_iteration = 0;
    md->UseAutoStep = false;
    md->time = 0.0;
    md->time_iteration = 0;

    md->SaveDirectory = (char *)m_alloc(1);
    md->SaveDirectory[0] = '\0';

    md->subd.grid = NULL;
    md->TotalNoOfParticles = 0;
    md->TotalNoOfMolecules = 0;
    md->ActiveGroup = FMD_GROUP_ALL;             /* all groups are active by default */
    md->ParticlesDistributed = false;
    md->ggrid = NULL;
    md->BoxSizeDetermined = false;
    md->PBCdetermined = false;
    md->PBC[0] = md->PBC[1] = md->PBC[2] = false;
    md->Is_MD_comm_root = false;
    md->Is_MD_process = false;
    md->EventHandler = NULL;
    md->timers = NULL;
    md->timers_num = 0;
    md->turies = NULL;
    md->turies_num = 0;
    md->active_ttm_turi = NULL;
    md->SaveConfigMode = FMD_SCM_XYZ_ATOMSNUM;
    md->cell_increment = 10;
    md->_OldNumberOfParticles = -1;
    md->_FileIndex = 0;
    md->KineticEnergyUpdated = false;
    md->random_seed_aux = (int)md;
    _fmd_potsys_init(md);
    create_mpi_types(md);
    _fmd_h5_ds_init(&md->h5_dataspaces);

    // this must be the last statement before return
    md->WallTimeOrigin = MPI_Wtime();

    return md;
}

void fmd_box_setSize(fmd_t *md, fmd_real_t sx, fmd_real_t sy, fmd_real_t sz)
{
    if (md->ggrid == NULL)
    {
        md->l[0] = sx;
        md->l[1] = sy;
        md->l[2] = sz;
        md->BoxSizeDetermined = true;
    }
}

fmd_real_t fmd_proc_getWallTime(fmd_t *md)
{
    return (MPI_Wtime() - md->WallTimeOrigin);
}

bool fmd_proc_hasSubdomain(fmd_t *md)
{
    return md->Is_MD_process;
}

bool fmd_proc_isRoot(fmd_t *md)
{
    return md->Is_MD_comm_root;
}

void fmd_proc_setNumThreads(fmd_t *md, int num)
{
    md->numthreads = num;
}

void fmd_io_setSaveDirectory(fmd_t *md, fmd_string_t directory)
{
    md->SaveDirectory = re_alloc(md->SaveDirectory, strlen(directory)+1);
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

void fmd_free(fmd_t *md)
{
    free(md->SaveDirectory);
    _fmd_subd_free(md);
    _fmd_potsys_free(md);
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

void fmd_setEventHandler(fmd_t *md, void *usp, fmd_EventHandler_t func)
{
    md->EventHandler = func;
    md->userobject = usp;
}
