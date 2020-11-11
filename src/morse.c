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
    fmd_ituple_t jc, kc;
    int d;
    particle_t *p1, *p2;
    fmd_real_t r2;
    fmd_rtuple_t rv;
    int ic0, ic1, ic2;
    fmd_real_t potEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    /* iterate over all cells */
    #pragma omp parallel for private(ic0,ic1,ic2,p1,d,kc,jc,p2,rv,r2) \
      shared(md,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
        for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
            for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
            {
                /* iterate over all particles in cell ic */
                for (int i1=0; i1 < md->SubDomain.grid[ic0][ic1][ic2].parts_num; i1++)
                {
                    p1 = &md->SubDomain.grid[ic0][ic1][ic2].parts[i1];

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
                                /* iterate over all particles in cell jc */
                                for (int i2=0; i2 < md->SubDomain.grid[jc[0]][jc[1]][jc[2]].parts_num; i2++)
                                {
                                    p2 = &md->SubDomain.grid[jc[0]][jc[1]][jc[2]].parts[i2];

                                    if (!(md->ActiveGroup == -1 || p2->core.GroupID == md->ActiveGroup))
                                        continue;

                                    if (p1 != p2)
                                    {
                                        unsigned atomkind2 = p2->core.atomkind;

                                        COMPUTE_rv_AND_r2;

                                        morse_t *morse = (morse_t *)pottable[atomkind1][atomkind2].data;

                                        if (r2 < morse->cutoff_sqr)
                                        {
                                            /* force, F = -(d/dr)U */
                                            fmd_real_t r = sqrt(r2);
                                            fmd_real_t inv_r = 1.0/r;
                                            fmd_real_t exp1 = exp( -morse->alpha * (r - morse->r0) );
                                            fmd_real_t exp2 = SQR(exp1);
                                            fmd_real_t factor = morse->alpha * morse->D0 * inv_r * (exp2 - exp1);
                                            for (d=0; d<3; d++)
                                                p1->F[d] += factor * rv[d];

                                            /* potential energy, U = D0 * ( exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)) ) */
                                            potEnergy += morse->D0 * (exp2 - 2.0 * exp1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (d=0; d<3; d++)
                        p1->F[d] *= 2;
                }
            }

    potEnergy *= 0.5;  /*correct fmd_real_t-counting*/
    MPI_Allreduce(&potEnergy, &md->TotalPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                               fmd_real_t D0, fmd_real_t alpha, fmd_real_t r0, fmd_real_t cutoff)
{
    morse_t *morse = (morse_t *)malloc(sizeof(morse_t));
    morse->D0 = D0;
    morse->alpha = alpha;
    morse->r0 = r0;
    morse->cutoff_sqr = SQR(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->cat = POT_MORSE;
    pot->data = morse;

    // add the pot to potlist
    md->potsys.potlist = fmd_list_prepend(md->potsys.potlist, pot);

    // apply the pot
    fmd_pot_apply(md, atomkind1, atomkind2, pot);

    return pot;
}
