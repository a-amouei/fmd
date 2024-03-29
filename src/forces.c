/*
  forces.c: This file is part of Free Molecular Dynamics

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

#include <tgmath.h>
#include "forces.h"
#include "fmd-private.h"
#include "eam.h"
#include "lj.h"
#include "morse.h"
#include "misc.h"
#include "md-ghost.h"
#include "list.h"
#include "molecule.h"
#include "general.h"
#include "turi.h"
#include "ttm.h"

static void compute_hybrid_pass1(fmd_t *md, fmd_real_t *FembSum_p)
{
    fmd_real_t Femb_sum = 0.0;
    potpair_t **pottable = md->potsys.pottable;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    /* iterate over all cells */

    #pragma omp parallel for shared(md,pottable,atomkinds) default(none) collapse(DIM) reduction(+:Femb_sum) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic0 = md->Subdomain.ic_start[0]; ic0 < md->Subdomain.ic_stop[0]; ic0++)
    for (int ic1 = md->Subdomain.ic_start[1]; ic1 < md->Subdomain.ic_stop[1]; ic1++)
    for (int ic2 = md->Subdomain.ic_start[2]; ic2 < md->Subdomain.ic_stop[2]; ic2++)
    {
        cell_t *c1 = &md->Subdomain.grid[ic0][ic1][ic2];

        /* iterate over all particles in cell ic */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            if (atomkinds[atomkind1].eam_element == NULL) continue;

            fmd_real_t rho_host = 0.0;

            /* iterate over neighbor cells of cell c1 */

            fmd_ituple_t jc, kc;

            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0);
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1);
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2);

                        cell_t *c2 = &ARRAY_ELEMENT(md->Subdomain.grid, jc);

                        // iterate over all items in cell c2

                        for (int i2=0; i2 < c2->parts_num; i2++)
                        {
                            if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                            if ( (c1 != c2) || (i1 != i2) )
                            {
                                unsigned atomkind2 = c2->atomkind[i2];

                                if (pottable[atomkind1][atomkind2].cat == POT_EAM_ALLOY)
                                    EAM_PAIR_UPDATE_rho_host(x1, atomkind1, kc, c1, i1, atomkind2,
                                                             c2, i2, rho_host, pottable);
                            }
                        }
                    }
                }
            }

            EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum(c1, i1, atomkind1, rho_host, atomkinds, Femb_sum);
        }
    }

    *FembSum_p = Femb_sum;
}

static void compute_hybrid_pass0(fmd_t *md, fmd_real_t FembSum)
{
    potpair_t **pottable = md->potsys.pottable;
    fmd_real_t PotEnergy = 0.0;

    /* iterate over all cells */

    #pragma omp parallel for shared(md,pottable) default(none) collapse(DIM) reduction(+:PotEnergy) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic0 = md->Subdomain.ic_start[0]; ic0 < md->Subdomain.ic_stop[0]; ic0++)
    for (int ic1 = md->Subdomain.ic_start[1]; ic1 < md->Subdomain.ic_stop[1]; ic1++)
    for (int ic2 = md->Subdomain.ic_start[2]; ic2 < md->Subdomain.ic_stop[2]; ic2++)
    {
        cell_t *c1 = &md->Subdomain.grid[ic0][ic1][ic2];

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            for (int d=0; d<DIM; d++)
                FRC(c1, i1, d) = 0.0;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            /* iterate over neighbor cells of cell c1 */

            fmd_ituple_t jc, kc;

            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0);
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1);
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2);

                        cell_t *c2 = &ARRAY_ELEMENT(md->Subdomain.grid, jc);

                        /* iterate over all particles in cell c2 */

                        for (int i2=0; i2 < c2->parts_num; i2++)
                        {
                            if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                            if ( (c1 != c2) || (i1 != i2) )
                            {
                                unsigned atomkind2 = c2->atomkind[i2];

                                switch (pottable[atomkind1][atomkind2].cat)
                                {
                                    case POT_EAM_ALLOY:
                                        EAM_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, kc, c1, i1, atomkind2,
                                                                            c2, i2, PotEnergy, pottable);
                                        break;

                                    case POT_LJ_6_12:
                                        LJ_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, kc, c1, i1, atomkind2,
                                                                           c2, i2, PotEnergy, pottable);
                                        break;

                                    case POT_MORSE:
                                        MORSE_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, kc, c1, i1, atomkind2,
                                                                              c2, i2, PotEnergy, pottable);
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    PotEnergy = 0.5 * PotEnergy + FembSum;
    MPI_Allreduce(&PotEnergy, &md->GroupPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

static void add_ttm_term_to_forces(fmd_t *md)
{
    fmd_ituple_t ic;
    cell_t *c;
    unsigned i;

    turi_t *t = md->active_ttm_turi;
    ttm_t *ttm = t->ttm;

    LOOP3D(ic, md->Subdomain.ic_start, md->Subdomain.ic_stop)
        for (c = &ARRAY_ELEMENT(md->Subdomain.grid, ic), i=0; i < c->parts_num; i++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c->GroupID[i] != md->ActiveGroup)
                continue;

            fmd_real_t mass = md->potsys.atomkinds[c->atomkind[i]].mass;

            if (ttm->dim == 1)
            {
                int itc1d = (int)(POS(c, i, DIM-1) / t->tcellh[DIM-1]) + t->itc_glob_to_loc[DIM-1];

                fmd_real_t xi = ttm->xi_1d[itc1d];
                fmd_real_t *vcm = ttm->vcm_1d[itc1d];

                FRC(c, i, DIM-1) += xi * mass * (VEL(c, i, DIM-1) - vcm[DIM-1]);
            }
            else
            {
                fmd_ituple_t itc;

                for (int d=0; d<DIM; d++)
                    itc[d] = (int)(POS(c, i, d) / t->tcellh[d]) + t->itc_glob_to_loc[d];

                fmd_real_t ximass = ARRAY_ELEMENT(ttm->xi, itc) * mass;
                fmd_real_t *vcm = ARRAY_ELEMENT(ttm->vcm, itc);

                for (int d=0; d<DIM; d++)
                    FRC(c, i, d) += ximass * (VEL(c, i, d) - vcm[d]);
            }
        }
}

void _fmd_dync_updateForces(fmd_t *md)
{
    if (md->potsys.potcats == NULL)  // just for one time
        _fmd_pot_prepareForForceComp(md);

    _fmd_ghostparticles_init(md);

    if (md->potsys.potcats_num == 1) // not hybrid mode
    {
        potcat_t potkind = *(potcat_t *)(md->potsys.potcats->data);

        switch (potkind)
        {
            case POT_LJ_6_12:
                fmd_computeLJ(md);
                break;

            case POT_MORSE:
                fmd_computeMorse(md);
                break;

            case POT_EAM_ALLOY: ;
                fmd_real_t FembSum;
                _fmd_computeEAM_pass1(md, &FembSum);
                _fmd_ghostparticles_update_Femb(md);
                _fmd_computeEAM_pass0(md, FembSum);
                break;
        }
    }
    else  // hybrid mode
    {
        fmd_real_t FembSum = 0.0;

        if (md->potsys.hybridpasses[1])
        {
            compute_hybrid_pass1(md, &FembSum);
            _fmd_ghostparticles_update_Femb(md);
        }

        if (md->potsys.hybridpasses[0])
            compute_hybrid_pass0(md, FembSum);
    }

    if (md->TotalNoOfMolecules > 0) fmd_dync_computeBondForce(md);   // TO-DO

    if (md->active_ttm_turi != NULL)
        if (_is_time_within_turi_start_stop_times(md, md->active_ttm_turi))
            add_ttm_term_to_forces(md);

    _fmd_ghostparticles_delete(md);
}
