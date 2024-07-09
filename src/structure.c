/*
  structure.c: This file is part of Free Molecular Dynamics

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

#include <time.h>
#include <tgmath.h>
#include <gsl/gsl_randist.h>
#include "fmd-private.h"
#include "misc.h"
#include "array.h"
#include "general.h"
#include "cell.h"

#define RESERVED_GROUP  -2

typedef enum {LATTICE_FCC, LATTICE_BCC, LATTICE_SC} lattice_t;

static void removeRemainingMomentum(fmd_t *md, int GroupID, fmd_real_t MomentumSum[], unsigned AtomsNum)
{
    cell_t *cell;
    int i;

    int nc = md->nc[0] * md->nc[1] * md->nc[2];

    for (int ic=0; ic < nc; ic++)
        for (cell = md->ggrid + ic, i = 0; i < cell->parts_num; i++)
            if (cell->GroupID[i] == RESERVED_GROUP)
            {
                cell->GroupID[i] = GroupID;

                fmd_real_t mass = md->potsys.atomkinds[cell->atomkind[i]].mass;

                for (int d=0; d<DIM; d++)
                    VEL(cell, i, d) -= MomentumSum[d] / (AtomsNum * mass);
            }
}

static void makeCuboid_mix(fmd_t *md, lattice_t lt, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp)
{
    const fmd_rtuple_t r_fcc[4] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                                   {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
    const fmd_rtuple_t r_bcc[2] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
    const fmd_rtuple_t r_sc[1] = {{0.0, 0.0, 0.0}};

    int points_per_cell;
    fmd_real_t *rp;

    switch (lt)
    {
        case LATTICE_SC:
            points_per_cell = 1;
            rp = (fmd_real_t *)r_sc;
            break;
        case LATTICE_FCC:
            points_per_cell = 4;
            rp = (fmd_real_t *)r_fcc;
            break;
        case LATTICE_BCC:
            points_per_cell = 2;
            rp = (fmd_real_t *)r_bcc;
            break;
    }

    fmd_ituple_t CrystalCell, ic;
    fmd_rtuple_t MomentumSum = {0.0, 0.0, 0.0};
    fmd_ituple_t dims = {dimx, dimy, dimz};
    fmd_rtuple_t r0 = {x, y, z};

    fmd_real_t *prps_cumult = (fmd_real_t *)m_alloc(md->potsys.atomkinds_num * sizeof(fmd_real_t));
    fmd_real_t prps_sum = 0.0;

    for (int i=0; i < md->potsys.atomkinds_num; i++)
        prps_cumult[i] = (prps_sum += ratio[i]);

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL) + md->random_seed_aux++);

    LOOP3D(CrystalCell, _fmd_ThreeZeros_int, dims)
        for (int i=0; i<points_per_cell; i++)
        {
            fmd_rtuple_t x;

            for (int d=0; d<DIM; d++)
            {
                x[d] = r0[d] + (CrystalCell[d] + .25 + rp[i*3+d]) * lp;
                ic[d] = (int)floor(x[d] / md->cellh[d]);
            }

            cell_t *c = md->ggrid + INDEX_FLAT(ic, md->nc);

            unsigned pi = _fmd_cell_new_particle(md, c);

            fmd_real_t rn = prps_sum * gsl_rng_uniform(rng);

            unsigned j;
            for (j=0; j < md->potsys.atomkinds_num; j++)
                if (rn < prps_cumult[j]) break;
            c->atomkind[pi] = j;

            fmd_real_t mass = md->potsys.atomkinds[j].mass;
            fmd_real_t StdDevVelocity = sqrt(K_BOLTZMANN * temp / mass);

            c->GroupID[pi] = RESERVED_GROUP;
            c->AtomID[pi] = md->TotalNoOfParticles++;

            for (int d=0; d<DIM; d++)
            {
                POS(c, pi, d) = x[d];
                VEL(c, pi, d) = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
                MomentumSum[d] += mass * VEL(c, pi, d);
            }
        }

    gsl_rng_free(rng);

    removeRemainingMomentum(md, GroupID, MomentumSum, points_per_cell * dimx * dimy * dimz);

    free(prps_cumult);
}

void fmd_matt_makeCuboidSC_mix(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp)
{
    if (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (GroupID == md->ActiveGroup || md->ActiveGroup == FMD_GROUP_ALL)
        md->KineticEnergyUpdated = false;

    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    makeCuboid_mix(md, LATTICE_SC, x, y, z, dimx, dimy, dimz, lp, ratio, GroupID, temp);
}

void fmd_matt_makeCuboidSC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, unsigned atomkind, int GroupID, fmd_real_t temp)
{
    if  (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (GroupID == md->ActiveGroup || md->ActiveGroup == FMD_GROUP_ALL)
        md->KineticEnergyUpdated = false;

    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    fmd_real_t *ratio = (fmd_real_t *)c_alloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    ratio[atomkind] = 1.0;

    makeCuboid_mix(md, LATTICE_SC, x, y, z, dimx, dimy, dimz, lp, ratio, GroupID, temp);

    free(ratio);
}

void fmd_matt_makeCuboidBCC_mix(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp)
{
    if  (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (GroupID == md->ActiveGroup || md->ActiveGroup == FMD_GROUP_ALL)
        md->KineticEnergyUpdated = false;

    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    makeCuboid_mix(md, LATTICE_BCC, x, y, z, dimx, dimy, dimz, lp, ratio, GroupID, temp);
}

void fmd_matt_makeCuboidBCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, unsigned atomkind, int GroupID, fmd_real_t temp)
{
    if  (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (GroupID == md->ActiveGroup || md->ActiveGroup == FMD_GROUP_ALL)
        md->KineticEnergyUpdated = false;

    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    fmd_real_t *ratio = (fmd_real_t *)c_alloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    ratio[atomkind] = 1.0;

    makeCuboid_mix(md, LATTICE_BCC, x, y, z, dimx, dimy, dimz, lp, ratio, GroupID, temp);

    free(ratio);
}

void fmd_matt_makeCuboidFCC_mix(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp)
{
    if (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (GroupID == md->ActiveGroup || md->ActiveGroup == FMD_GROUP_ALL)
        md->KineticEnergyUpdated = false;

    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    makeCuboid_mix(md, LATTICE_FCC, x, y, z, dimx, dimy, dimz, lp, ratio, GroupID, temp);
}

void fmd_matt_makeCuboidFCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, unsigned atomkind, int GroupID, fmd_real_t temp)
{
    if (md->ggrid == NULL) _fmd_createGlobalGrid(md);

    if (GroupID == md->ActiveGroup || md->ActiveGroup == FMD_GROUP_ALL)
        md->KineticEnergyUpdated = false;

    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    fmd_real_t *ratio = (fmd_real_t *)c_alloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    ratio[atomkind] = 1.0;

    makeCuboid_mix(md, LATTICE_FCC, x, y, z, dimx, dimy, dimz, lp, ratio, GroupID, temp);

    free(ratio);
}
