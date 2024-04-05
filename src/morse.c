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

#include <tgmath.h>
#include "morse.h"
#include "fmd-private.h"
#include "misc.h"
#include "forces.h"
#include "list.h"
#include "potential.h"
#include "general.h"
#include "cell.h"

void fmd_computeMorse(fmd_t *md)
{
    fmd_real_t PotEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    /* iterate over all cells */

    #pragma omp parallel for shared(md,pottable) default(none) reduction(+:PotEnergy) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic=0; ic < md->subd.nc; ic++)
    {
        cell_t *c1 = md->subd.grid + ic;

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);
            fmd_real_t *F1 = &FRC(c1, i1, 0);

            for (int d=0; d < DIM; d++)
                F1[d] = 0.0;

            /* iterate over neighbor cells of cell c1 */

            for (int jc=0; jc < CNEIGHBS_NUM; jc++)
            {
                cell_t *c2 = c1->cneighbs[jc];

                /* iterate over all particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    if ( (c1 != c2) || (i1 != i2) )
                    {
                        fmd_real_t r2;
                        fmd_rtuple_t rv;

                        unsigned atomkind2 = c2->atomkind[i2];
                        fmd_real_t *x2 = &POS(c2, i2, 0);

                        COMPUTE_rv_AND_r2(x1, x2, rv, r2);

                        morse_t *morse = (morse_t *)pottable[atomkind1][atomkind2].data;

                        if (r2 < morse->cutoff_sqr)
                        {
                            /* force, F = -(d/dr)U */

                            fmd_real_t r = sqrt(r2);
                            fmd_real_t inv_r = 1.0/r;
                            fmd_real_t exp1 = exp( -morse->alpha * (r - morse->r0) );
                            fmd_real_t exp2 = sqrr(exp1);
                            fmd_real_t factor = morse->alpha * morse->D0 * inv_r * (exp2 - exp1);

                            for (int d=0; d<DIM; d++)
                                F1[d] += factor * rv[d];

                            /* potential energy, U = D0 * ( exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)) ) */
                            PotEnergy += morse->D0 * (exp2 - 2.0 * exp1);
                        }
                    }
                }
            }

            for (int d=0; d<DIM; d++)
                F1[d] *= 2;
        }
    }

    PotEnergy *= 0.5;

    MPI_Allreduce(&PotEnergy, &md->GroupPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                               fmd_real_t D0, fmd_real_t alpha, fmd_real_t r0, fmd_real_t cutoff)
{
    morse_t *morse = (morse_t *)m_alloc(sizeof(morse_t));
    morse->D0 = D0;
    morse->alpha = alpha;
    morse->r0 = r0;
    morse->cutoff_sqr = sqrr(cutoff);

    fmd_pot_t *pot = (fmd_pot_t *)m_alloc(sizeof(fmd_pot_t));
    pot->cat = POT_MORSE;
    pot->data = morse;

    // add the pot to potlist
    md->potsys.potlist = _fmd_list_prepend(md->potsys.potlist, pot);

    // apply the pot
    fmd_pot_apply(md, atomkind1, atomkind2, pot);

    return pot;
}
