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
    fmd_ituple_t jc, kc;
    int d;
    fmd_real_t r2;
    fmd_rtuple_t rv;
    int ic0, ic1, ic2;
    fmd_real_t potEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    /* iterate over all cells(lists) */
    #pragma omp parallel for private(ic0,ic1,ic2,d,kc,jc,rv,r2) \
      shared(md,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
        for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
            for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
            {
                int i1, i2;
                cell_t *cell1, *cell2;

                /* iterate over all items in cell ic */
                for (cell1 = &md->SubDomain.grid[ic0][ic1][ic2], i1=0; i1 < cell1->parts_num; i1++)
                {
                    particle_t *p1 = &cell1->parts[i1];

                    if (!(md->ActiveGroup == -1 || p1->core.GroupID == md->ActiveGroup))
                        continue;

                    unsigned atomkind1 = p1->core.atomkind;

                    for (d=0; d<3; d++)
                        p1->F[d] = 0.0;

                    /* iterate over neighbor cells of cell ic */
                    for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
                    {
                        SET_jc_IN_DIRECTION(0)
                        for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                        {
                            SET_jc_IN_DIRECTION(1)
                            for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                            {
                                SET_jc_IN_DIRECTION(2)

                                /* iterate over all items in cell jc */
                                for (cell2 = &md->SubDomain.grid[jc[0]][jc[1]][jc[2]], i2=0; i2 < cell2->parts_num; i2++)
                                {
                                    particle_t *p2 = &cell2->parts[i2];

                                    if (!(md->ActiveGroup == -1 || p2->core.GroupID == md->ActiveGroup))
                                        continue;

                                    if (p1->core.molkind!=0 && p1->core.MolID==p2->core.MolID) continue;// TO-DO

                                    if (p1 != p2)
                                    {
                                        unsigned atomkind2 = p2->core.atomkind;

                                        COMPUTE_rv_AND_r2;

                                        LJ_6_12_t *lj = (LJ_6_12_t *)pottable[atomkind1][atomkind2].data;

                                        if (r2 < lj->cutoff_sqr)
                                        {
                                            fmd_real_t inv_r2, inv_rs2, inv_rs6, inv_rs12;

                                            /* force, F = -(d/dr)U */
                                            inv_r2 = 1.0/r2;
                                            inv_rs2 = SQR(lj->sig) * inv_r2;
                                            inv_rs6 = inv_rs2 * inv_rs2 * inv_rs2;
                                            inv_rs12 = SQR(inv_rs6);
                                            fmd_real_t factor = lj->eps * inv_r2 * (inv_rs12 - 0.5*inv_rs6);
                                            for (d=0; d<3; d++)
                                                p1->F[d] += rv[d] * factor;

                                            /* potential energy, U = 4*eps*( (sig/r)^12 - (sig/r)^6 ) */
                                            potEnergy += lj->eps * (inv_rs12 - inv_rs6);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    for (d=0; d<3; d++)
                        p1->F[d] *= 48.0;
                }
            }

    potEnergy *= 2.0;
    MPI_Allreduce(&potEnergy, &md->TotalPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                            fmd_real_t sigma, fmd_real_t epsilon, fmd_real_t cutoff)
{
    LJ_6_12_t *lj = (LJ_6_12_t *)malloc(sizeof(LJ_6_12_t));
    lj->sig = sigma;
    lj->eps = epsilon;
    lj->cutoff_sqr = SQR(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->cat = POT_LJ_6_12;
    pot->data = lj;

    /* add the pot to potlist */
    md->potsys.potlist = fmd_list_prepend(md->potsys.potlist, pot);

    /* apply the pot */
    fmd_pot_apply(md, atomkind1, atomkind2, pot);

    return pot;
}
