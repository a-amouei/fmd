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

static void removeRemainingMomentum(fmd_t *md, int GroupID, double *MomentumSum, unsigned AtomsNum)
{
    TParticleListItem *item_p;
    int ic[3];

    ITERATE(ic, fmd_ThreeZeros, md->nc)
        for (item_p=md->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            if (item_p->P.GroupID == RESERVED_GROUP)
            {
                item_p->P.GroupID = GroupID;
                double mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
                for (int d=0; d<3; d++)
                    item_p->P.v[d] -= MomentumSum[d] / (AtomsNum * mass);
            }
        }
}

void fmd_matt_makeCuboidFCC_alloy(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, double *proportions, int GroupID)
{
    if (md->Is_MD_comm_root)
    {
        double mass, StdDevVelocity;
        int CrystalCell[3], ic[3];
        double r_fcc[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                              {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
        double MomentumSum[3] = {0.0, 0.0, 0.0};
        TParticleListItem *item_p;
        int dims[3] = {dimx, dimy, dimz};
        double r0[3] = {x, y, z};
        int i, d;

        double *prps_cumult = (double *)malloc(md->potsys.atomkinds_num * sizeof(double));
        double prps_sum = 0.0;
        for (i=0; i < md->potsys.atomkinds_num; i++)
            prps_cumult[i] = (prps_sum += proportions[i]);

        gsl_rng *rng;
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, time(NULL));

        ITERATE(CrystalCell, fmd_ThreeZeros, dims)
            for (i=0; i<4; i++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));

                double rn = prps_sum * gsl_rng_uniform(rng);
                int j;
                for (j=0; j < md->potsys.atomkinds_num; j++)
                    if (rn < prps_cumult[j]) break;
                item_p->P.atomkind = j;

                mass = md->potsys.atomkinds[j].mass;
                StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

                item_p->P.GroupID = RESERVED_GROUP;
                item_p->P.AtomID = md->TotalNoOfParticles++;

                for (d=0; d<3; d++)
                {
                    item_p->P.x[d] = r0[d] + (CrystalCell[d] + .25 + r_fcc[i][d]) * LatticeParameter;
                    item_p->P.v[d] = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
                    MomentumSum[d] += mass * item_p->P.v[d];
                    ic[d] = (int)floor(item_p->P.x[d] / md->cellh[d]);
                }
                fmd_insertInList(&md->global_grid[ic[0]][ic[1]][ic[2]], item_p);
            }

        gsl_rng_free(rng);

        removeRemainingMomentum(md, GroupID, MomentumSum, 4 * dimx * dimy * dimz);

        free(prps_cumult);
    }
}

void fmd_matt_makeCuboidFCC(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    double *proportions = (double *)calloc(md->potsys.atomkinds_num, sizeof(double));
    proportions[atomkind] = 1.0;

    fmd_matt_makeCuboidFCC_alloy(md, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_scatterMolecule(fmd_t *md, unsigned molkind, double xa,
  double ya, double za, double xb, double yb, double zb, unsigned num,
  int GroupID)
{
    if (!md->Is_MD_comm_root) return;

    double x1[3] = {xa, ya, za};
    double x2[3] = {xb, yb, zb};

    double X[3];
    for (int d=0; d<3; d++)
        X[d] = x2[d] - x1[d];

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    molkind_t *mk = &md->potsys.molkinds[molkind];
    double **coords = (double **)fmd_array_neat2d_create(mk->atoms_num, 3, sizeof(double));
    double MomentumSum[3] = {0., 0., 0.};

    for (unsigned i=0; i<num; i++) // iteration on number of molecules
    {
        unsigned j, trycount=0;

        // prepare coordinates for this molecule
        do
        {
            double xo[3];

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

        TParticleListItem *item_p;
        int ic[3];

        // create the atoms in memory and initialize them
        for (j=0; j < mk->atoms_num; j++)
        {
            item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
            item_p->P.AtomID_local = j;
            item_p->P.atomkind = mk->atoms[j].atomkind;
            item_p->P.MolID = md->TotalNoOfMolecules;
            item_p->P.molkind = molkind;
            item_p->P.AtomID = md->TotalNoOfParticles++;
            item_p->P.GroupID = RESERVED_GROUP;

            double mass = md->potsys.atomkinds[item_p->P.atomkind].mass;
            double StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

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

    fmd_array_neat2d_free((void **)coords);
}
