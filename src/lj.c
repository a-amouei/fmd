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
#include "potential.h"
#include "base.h"
#include "list.h"
#include "forces.h"

void fmd_computeLJ(fmd_t *md)
{
    int jc[3], kc[3];
    int d;
    ParticleListItem_t *item1_p, *item2_p;
    double r2, rv[3];
    int ic0, ic1, ic2;
    double potEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,kc,jc,item2_p,rv,r2) \
      shared(md,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
        for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
            for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
            {
                // iterate over all items in cell ic
                for (item1_p = md->SubDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
                {
                    if (!(md->ActiveGroup == -1 || item1_p->P.GroupID == md->ActiveGroup))
                        continue;

                    unsigned atomkind1 = item1_p->P.atomkind;

                    for (d=0; d<3; d++)
                        item1_p->F[d] = 0.0;

                    // iterate over neighbor cells of cell ic
                    for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
                    {
                        SET_jc_IN_DIRECTION(0)
                        for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                        {
                            SET_jc_IN_DIRECTION(1)
                            for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                            {
                                SET_jc_IN_DIRECTION(2)
                                // iterate over all items in cell jc
                                for (item2_p = md->SubDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                                {
                                    if (!(md->ActiveGroup == -1 || item2_p->P.GroupID == md->ActiveGroup))
                                        continue;

                                    if (item1_p->P.molkind!=0 && item1_p->P.MolID==item2_p->P.MolID) continue;// TO-DO

                                    if (item1_p != item2_p)
                                    {
                                        unsigned atomkind2 = item2_p->P.atomkind;

                                        COMPUTE_rv_AND_r2;

                                        LJ_6_12_t *lj = (LJ_6_12_t *)pottable[atomkind1][atomkind2].data;

                                        if (r2 < lj->cutoff_sqr)
                                        {
                                            double inv_r2, inv_rs2, inv_rs6, inv_rs12;

                                            // force, F = -(d/dr)U
                                            inv_r2 = 1.0/r2;
                                            inv_rs2 = SQR(lj->sig) * inv_r2;
                                            inv_rs6 = inv_rs2 * inv_rs2 * inv_rs2;
                                            inv_rs12 = SQR(inv_rs6);
                                            double factor = lj->eps * inv_r2 * (inv_rs12 - 0.5*inv_rs6);
                                            for (d=0; d<3; d++)
                                                item1_p->F[d] += rv[d] * factor;

                                            // potential energy, U = 4*eps*( (sig/r)^12 - (sig/r)^6 )
                                            potEnergy += lj->eps * (inv_rs12 - inv_rs6);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    for (d=0; d<3; d++)
                        item1_p->F[d] *= 48.0;
                }
            }

    potEnergy *= 2.0;
    MPI_Allreduce(&potEnergy, &md->TotalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                            double sigma, double epsilon, double cutoff)
{
    LJ_6_12_t *lj = (LJ_6_12_t *)malloc(sizeof(LJ_6_12_t));
    lj->sig = sigma;
    lj->eps = epsilon;
    lj->cutoff_sqr = SQR(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->cat = POT_LJ_6_12;
    pot->data = lj;

    // add the pot to potlist
    md->potsys.potlist = fmd_list_prepend(md->potsys.potlist, pot);

    // apply the pot
    fmd_pot_apply(md, atomkind1, atomkind2, pot);

    return pot;
}
