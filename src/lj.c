/*
  lj.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

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

#include "lj.h"
#include "fmd-private.h"
#include "potential.h"
#include "misc.h"
#include "list.h"
#include "forces.h"
#include "general.h"

void fmd_computeLJ(fmd_t *md)
{
    fmd_real_t PotEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    /* iterate over all cells (lists) */

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

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);
            fmd_real_t *F1 = &FRC(c1, i1, 0);

            for (int d=0; d<DIM; d++)
                F1[d] = 0.0;

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

                            /*
                            if (c1->molkind != NULL)
                                if (c1->molkind[i1] != 0 && c1->MolID[i1] == c2->MolID[i2]) continue;  // TO-DO
                            */

                            if ( (c1 != c2) || (i1 != i2) )
                            {
                                fmd_real_t r2;
                                fmd_rtuple_t rv;

                                unsigned atomkind2 = c2->atomkind[i2];
                                fmd_real_t *x2 = &POS(c2, i2, 0);

                                COMPUTE_rv_AND_r2(x1, x2, kc, rv, r2);

                                LJ_6_12_t *lj = (LJ_6_12_t *)pottable[atomkind1][atomkind2].data;

                                if (r2 < lj->cutoff_sqr)
                                {
                                    fmd_real_t inv_r2, inv_rs2, inv_rs6, inv_rs12;

                                    /* force, F = -(d/dr)U */
                                    inv_r2 = 1.0/r2;
                                    inv_rs2 = sqrr(lj->sig) * inv_r2;
                                    inv_rs6 = inv_rs2 * inv_rs2 * inv_rs2;
                                    inv_rs12 = sqrr(inv_rs6);
                                    fmd_real_t factor = lj->eps * inv_r2 * (inv_rs12 - 0.5*inv_rs6);

                                    for (int d=0; d<DIM; d++)
                                        F1[d] += rv[d] * factor;

                                    /* potential energy, U = 4*eps*( (sig/r)^12 - (sig/r)^6 ) */
                                    PotEnergy += lj->eps * (inv_rs12 - inv_rs6);
                                }
                            }
                        }
                    }
                }
            }

            for (int d=0; d<DIM; d++)
                F1[d] *= 48.0;
        }
    }

    PotEnergy *= 2.0;
    MPI_Allreduce(&PotEnergy, &md->GroupPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                            fmd_real_t sigma, fmd_real_t epsilon, fmd_real_t cutoff)
{
    LJ_6_12_t *lj = (LJ_6_12_t *)m_alloc(sizeof(LJ_6_12_t));
    lj->sig = sigma;
    lj->eps = epsilon;
    lj->cutoff_sqr = sqrr(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)m_alloc(sizeof(fmd_pot_t));
    pot->cat = POT_LJ_6_12;
    pot->data = lj;

    /* add the pot to potlist */
    md->potsys.potlist = _fmd_list_prepend(md->potsys.potlist, pot);

    /* apply the pot */
    fmd_pot_apply(md, atomkind1, atomkind2, pot);

    return pot;
}
