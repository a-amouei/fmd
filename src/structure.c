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

#include "base.h"
#include "molecule.h"
#include "array.h"

#define RESERVED_GROUP  -2

typedef enum {LATTICE_FCC, LATTICE_BCC, LATTICE_SC} lattice_t;

static void removeRemainingMomentum(fmd_t *md, int GroupID, fmd_real_t *MomentumSum, unsigned AtomsNum)
{
    ParticleListItem_t *item_p;
    fmd_ituple_t ic;

    ITERATE(ic, fmd_ThreeZeros, md->nc)
        for (item_p=md->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            if (item_p->P.GroupID == RESERVED_GROUP)
            {
                item_p->P.GroupID = GroupID;
                fmd_real_t mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
                for (int d=0; d<3; d++)
                    item_p->P.v[d] -= MomentumSum[d] / (AtomsNum * mass);
            }
        }
}

static void fmd_matt_makeCuboid_alloy(fmd_t *md, lattice_t lt, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    fmd_real_t mass, StdDevVelocity;
    fmd_ituple_t CrystalCell, ic;

    const fmd_rtuple_t r_fcc[4] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                                   {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
    const fmd_rtuple_t r_bcc[2] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
    const fmd_rtuple_t r_sc[1] = {{0.0, 0.0, 0.0}};

    fmd_rtuple_t MomentumSum = {0.0, 0.0, 0.0};
    ParticleListItem_t *item_p;
    fmd_ituple_t dims = {dimx, dimy, dimz};
    fmd_rtuple_t r0 = {x, y, z};
    int i, d;

    fmd_real_t *prps_cumult = (fmd_real_t *)malloc(md->potsys.atomkinds_num * sizeof(fmd_real_t));
    fmd_real_t prps_sum = 0.0;
    for (i=0; i < md->potsys.atomkinds_num; i++)
        prps_cumult[i] = (prps_sum += proportions[i]);

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

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

    ITERATE(CrystalCell, fmd_ThreeZeros, dims)
        for (i=0; i<points_per_cell; i++)
        {
            item_p = (ParticleListItem_t *)malloc(sizeof(ParticleListItem_t));

            fmd_real_t rn = prps_sum * gsl_rng_uniform(rng);
            int j;
            for (j=0; j < md->potsys.atomkinds_num; j++)
                if (rn < prps_cumult[j]) break;
            item_p->P.atomkind = j;

            item_p->P.molkind = 0;

            mass = md->potsys.atomkinds[j].mass;
            StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

            item_p->P.GroupID = RESERVED_GROUP;
            item_p->P.AtomID = md->TotalNoOfParticles++;

            for (d=0; d<3; d++)
            {
                item_p->P.x[d] = r0[d] + (CrystalCell[d] + .25 + rp[i*3+d]) * LatticeParameter;
                item_p->P.v[d] = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
                MomentumSum[d] += mass * item_p->P.v[d];
                ic[d] = (int)floor(item_p->P.x[d] / md->cellh[d]);
            }
            fmd_insertInList(&md->global_grid[ic[0]][ic[1]][ic[2]], item_p);
        }

    gsl_rng_free(rng);

    removeRemainingMomentum(md, GroupID, MomentumSum, points_per_cell * dimx * dimy * dimz);

    free(prps_cumult);
}

void fmd_matt_makeCuboidSC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_matt_makeCuboid_alloy(md, LATTICE_SC, x, y, z, dimx, dimy, dimz,
      LatticeParameter, proportions, GroupID);
}

void fmd_matt_makeCuboidSC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_real_t *proportions = (fmd_real_t *)calloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    proportions[atomkind] = 1.0;

    fmd_matt_makeCuboid_alloy(md, LATTICE_SC, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_makeCuboidBCC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_matt_makeCuboid_alloy(md, LATTICE_BCC, x, y, z, dimx, dimy, dimz,
      LatticeParameter, proportions, GroupID);
}

void fmd_matt_makeCuboidBCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_real_t *proportions = (fmd_real_t *)calloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    proportions[atomkind] = 1.0;

    fmd_matt_makeCuboid_alloy(md, LATTICE_BCC, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_makeCuboidFCC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_matt_makeCuboid_alloy(md, LATTICE_FCC, x, y, z, dimx, dimy, dimz,
      LatticeParameter, proportions, GroupID);
}

void fmd_matt_makeCuboidFCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_real_t *proportions = (fmd_real_t *)calloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    proportions[atomkind] = 1.0;

    fmd_matt_makeCuboid_alloy(md, LATTICE_FCC, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_scatterMolecule(fmd_t *md, unsigned molkind, fmd_real_t xa,
  fmd_real_t ya, fmd_real_t za, fmd_real_t xb, fmd_real_t yb, fmd_real_t zb, unsigned num,
  int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    fmd_rtuple_t x1 = {xa, ya, za};
    fmd_rtuple_t x2 = {xb, yb, zb};

    fmd_rtuple_t X;
    for (int d=0; d<3; d++)
        X[d] = x2[d] - x1[d];

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    molkind_t *mk = &md->potsys.molkinds[molkind];
    fmd_real_t **coords = (fmd_real_t **)_fmd_array_neat2d_create(mk->atoms_num, 3, sizeof(fmd_real_t));
    fmd_rtuple_t MomentumSum = {0., 0., 0.};

    for (unsigned i=0; i<num; i++) // iteration on number of molecules
    {
        unsigned j, trycount=0;

        // prepare coordinates for this molecule
        do
        {
            fmd_rtuple_t xo;

            for (int d=0; d<3; d++)  // choose center position
                xo[d] = x1[d] + X[d] * gsl_rng_uniform(rng);

            for (j=0; j < mk->atoms_num; j++)  // iteration on number of atoms in a molecule
            {
                fmd_bool_t broke = 0;

                for (int d=0; d<3; d++)
                {
                    coords[j][d] = xo[d] + mk->atoms[j].position[d];
                    if (coords[j][d] <= x1[d])
                    {
                        broke = 1;
                        break;
                    }
                    else if (coords[j][d] >= x2[d])
                    {
                        broke = 1;
                        break;
                    }
                }

                if (broke) break;
            }
            trycount++;
        } while (j != mk->atoms_num && trycount < 1000);
        // TO-DO: error handling
        assert(trycount < 1000);

        ParticleListItem_t *item_p;
        fmd_ituple_t ic;

        // create the atoms in memory and initialize them
        for (j=0; j < mk->atoms_num; j++)
        {
            item_p = (ParticleListItem_t *)malloc(sizeof(ParticleListItem_t));
            item_p->P.AtomID_local = j;
            item_p->P.atomkind = mk->atoms[j].atomkind;
            item_p->P.MolID = md->TotalNoOfMolecules;
            item_p->P.molkind = molkind;
            item_p->P.AtomID = md->TotalNoOfParticles++;
            item_p->P.GroupID = RESERVED_GROUP;

            fmd_real_t mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
            fmd_real_t StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

            for (int d=0; d<3; d++)
            {
                item_p->P.x[d] = coords[j][d];
                item_p->P.v[d] = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
                MomentumSum[d] += mass * item_p->P.v[d];
                ic[d] = (int)floor(item_p->P.x[d] / md->cellh[d]);
            }
            fmd_insertInList(&md->global_grid[ic[0]][ic[1]][ic[2]], item_p);
        }

        md->TotalNoOfMolecules++;
    }

    removeRemainingMomentum(md, GroupID, MomentumSum, num * mk->atoms_num);

    _fmd_array_neat2d_free((void **)coords);
}
