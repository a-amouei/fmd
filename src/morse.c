/*
  morse.c: This file is part of Free Molecular Dynamics

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

#include "morse.h"
#include "base.h"
#include "forces.h"
#include "list.h"
#include "potential.h"

void fmd_computeMorse(fmd_t *md)
{
    int jc[3], kc[3];
    int d;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    int ic0, ic1, ic2;
    double potEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,kc,jc,item2_p,rv,r2) \
      shared(md,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
    for (ic0 = md->subDomain.ic_start[0]; ic0 < md->subDomain.ic_stop[0]; ic0++)
        for (ic1 = md->subDomain.ic_start[1]; ic1 < md->subDomain.ic_stop[1]; ic1++)
            for (ic2 = md->subDomain.ic_start[2]; ic2 < md->subDomain.ic_stop[2]; ic2++)
            {
                // iterate over all items in cell ic
                for (item1_p = md->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
                {
                    if (!(md->activeGroup == -1 || item1_p->P.groupID == md->activeGroup))
                        continue;

                    unsigned atomkind1 = item1_p->P.elementID;

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
                                for (item2_p = md->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                                {
                                    if (!(md->activeGroup == -1 || item2_p->P.groupID == md->activeGroup))
                                        continue;

                                    if (item1_p != item2_p)
                                    {
                                        unsigned atomkind2 = item2_p->P.elementID;

                                        COMPUTE_rv_AND_r2;

                                        morse_t *morse = (morse_t *)pottable[atomkind1][atomkind2].data;

                                        if (r2 < morse->cutoff_sqr)
                                        {
                                            // force, F = -(d/dr)U
                                            double r = sqrt(r2);
                                            double inv_r = 1.0/r;
                                            double exp1 = exp( -morse->alpha * (r - morse->r0) );
                                            double exp2 = SQR(exp1);
                                            double factor = morse->alpha * morse->D0 * inv_r * (exp2 - exp1);
                                            for (d=0; d<3; d++)
                                                item1_p->F[d] += factor * rv[d];

                                            // potential energy, U = D0 * ( exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)) )
                                            potEnergy += morse->D0 * (exp2 - 2.0 * exp1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (d=0; d<3; d++)
                        item1_p->F[d] *= 2;
                }
            }

    potEnergy *= 0.5;  /*correct double-counting*/
    MPI_Allreduce(&potEnergy, &md->totalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                               double D0, double alpha, double r0, double cutoff)
{
    morse_t *morse = (morse_t *)malloc(sizeof(morse_t));
    morse->D0 = D0;
    morse->alpha = alpha;
    morse->r0 = r0;
    morse->cutoff_sqr = SQR(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->kind = POTKIND_MORSE;
    pot->data = morse;

    // add the pot to potlist
    md->potsys.potlist = fmd_list_prepend(md->potsys.potlist, pot);

    // apply the pot
    fmd_pot_apply(md, atomkind1, atomkind2, pot);

    return pot;
}
