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
#include "molecule.h"
#include "array.h"
#include "general.h"
#include "cell.h"

#define RESERVED_GROUP  -2

typedef enum {LATTICE_FCC, LATTICE_BCC, LATTICE_SC} lattice_t;

static void removeRemainingMomentum(fmd_t *md, int GroupID, fmd_real_t *MomentumSum, unsigned AtomsNum)
{
    fmd_ituple_t ic;
    cell_t *cell;
    int i;

    LOOP3D(ic, _fmd_ThreeZeros_int, md->nc)
        for (cell = &ARRAY_ELEMENT(md->global_grid, ic), i = 0; i < cell->parts_num; i++)
            if (cell->GroupID[i] == RESERVED_GROUP)
            {
                cell->GroupID[i] = GroupID;

                fmd_real_t mass = md->potsys.atomkinds[cell->atomkind[i]].mass;

                for (int d=0; d<3; d++)
                    VEL(cell, i, d) -= MomentumSum[d] / (AtomsNum * mass);
            }
}

static void makeCuboid_alloy(fmd_t *md, lattice_t lt, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
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
        prps_cumult[i] = (prps_sum += proportions[i]);

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    LOOP3D(CrystalCell, _fmd_ThreeZeros_int, dims)
        for (int i=0; i<points_per_cell; i++)
        {
            fmd_rtuple_t x;

            for (int d=0; d<3; d++)
            {
                x[d] = r0[d] + (CrystalCell[d] + .25 + rp[i*3+d]) * LatticeParameter;
                ic[d] = (int)floor(x[d] / md->cellh[d]);
            }

            cell_t *c = &ARRAY_ELEMENT(md->global_grid, ic);

            unsigned pi = _fmd_cell_new_particle(md, c);

            fmd_real_t rn = prps_sum * gsl_rng_uniform(rng);

            unsigned j;
            for (j=0; j < md->potsys.atomkinds_num; j++)
                if (rn < prps_cumult[j]) break;
            c->atomkind[pi] = j;

            if (c->molkind != NULL) c->molkind[pi] = 0;

            fmd_real_t mass = md->potsys.atomkinds[j].mass;
            fmd_real_t StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

            c->GroupID[pi] = RESERVED_GROUP;
            c->AtomID[pi] = md->TotalNoOfParticles++;

            for (int d=0; d<3; d++)
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

void fmd_matt_makeCuboidSC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    makeCuboid_alloy(md, LATTICE_SC, x, y, z, dimx, dimy, dimz,
      LatticeParameter, proportions, GroupID);
}

void fmd_matt_makeCuboidSC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    fmd_real_t *proportions = (fmd_real_t *)c_alloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    proportions[atomkind] = 1.0;

    makeCuboid_alloy(md, LATTICE_SC, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_makeCuboidBCC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    makeCuboid_alloy(md, LATTICE_BCC, x, y, z, dimx, dimy, dimz,
      LatticeParameter, proportions, GroupID);
}

void fmd_matt_makeCuboidBCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    fmd_real_t *proportions = (fmd_real_t *)c_alloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    proportions[atomkind] = 1.0;

    makeCuboid_alloy(md, LATTICE_BCC, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_makeCuboidFCC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    makeCuboid_alloy(md, LATTICE_FCC, x, y, z, dimx, dimy, dimz,
      LatticeParameter, proportions, GroupID);
}

void fmd_matt_makeCuboidFCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

    fmd_real_t *proportions = (fmd_real_t *)c_alloc(md->potsys.atomkinds_num, sizeof(fmd_real_t));
    proportions[atomkind] = 1.0;

    makeCuboid_alloy(md, LATTICE_FCC, x, y, z, dimx, dimy, dimz, LatticeParameter, proportions, GroupID);

    free(proportions);
}

void fmd_matt_scatterMolecule(fmd_t *md, fmd_handle_t molkind, fmd_real_t xa,
  fmd_real_t ya, fmd_real_t za, fmd_real_t xb, fmd_real_t yb, fmd_real_t zb, unsigned num,
  int GroupID)
{
    if (!md->GlobalGridExists) _fmd_createGlobalGrid(md);
    if (!md->Is_MD_comm_root) return;

    assert(GroupID >= 0); /* TO-DO: handle error */

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

    for (unsigned i=0; i<num; i++) /* iteration on number of molecules */
    {
        unsigned j, trycount=0;

        /* prepare coordinates for this molecule */
        do
        {
            fmd_rtuple_t xo;

            for (int d=0; d<3; d++)  /* choose center position */
                xo[d] = x1[d] + X[d] * gsl_rng_uniform(rng);

            for (j=0; j < mk->atoms_num; j++)  /* iteration on number of atoms in a molecule */
            {
                bool broke = false;

                for (int d=0; d<3; d++)
                {
                    coords[j][d] = xo[d] + mk->atoms[j].position[d];
                    if (coords[j][d] <= x1[d])
                    {
                        broke = true;
                        break;
                    }
                    else if (coords[j][d] >= x2[d])
                    {
                        broke = true;
                        break;
                    }
                }

                if (broke) break;
            }
            trycount++;
        } while (j != mk->atoms_num && trycount < 1000);
        /* TO-DO: error handling */
        assert(trycount < 1000);

        fmd_ituple_t ic;

        /* create the atoms in memory and initialize them */
        for (j=0; j < mk->atoms_num; j++)
        {
            fmd_rtuple_t x;

            for (int d=0; d<3; d++)
            {
                x[d] = coords[j][d];
                ic[d] = (int)floor(x[d] / md->cellh[d]);
            }

            cell_t *c = &ARRAY_ELEMENT(md->global_grid, ic);
            unsigned pi = _fmd_cell_new_particle(md, c);

            c->AtomIDlocal[pi] = j;
            c->atomkind[pi] = mk->atoms[j].atomkind;
            c->MolID[pi] = md->TotalNoOfMolecules;
            c->molkind[pi] = molkind;
            c->AtomID[pi] = md->TotalNoOfParticles++;
            c->GroupID[pi] = RESERVED_GROUP;

            fmd_real_t mass = md->potsys.atomkinds[c->atomkind[pi]].mass;
            fmd_real_t StdDevVelocity = sqrt(K_BOLTZMANN * md->DesiredTemperature / mass);

            for (int d=0; d<3; d++)
            {
                POS(c, pi, d) = x[d];
                VEL(c, pi, d) = gsl_ran_gaussian_ziggurat(rng, StdDevVelocity);
                MomentumSum[d] += mass * VEL(c, pi, d);
            }
        }

        md->TotalNoOfMolecules++;
    }

    removeRemainingMomentum(md, GroupID, MomentumSum, num * mk->atoms_num);

    _fmd_array_neat2d_free((void **)coords);
}
