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

void _fmd_computeLJ(fmd_t *md)
{
    fmd_real_t PotEnergy = 0.0;
    potpair_t **pottable = md->potsys.pottable;

    _fmd_clean_forces(md);

    /* iterate over all cells (lists) */

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

            /* iterate over neighbor cells of cell c1 */

            for (int jc = 0; jc < c1->cnb1len; jc++)
            {
                cell_t *c2 = c1->cnb1[jc];

                /* iterate over all particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    /*
                    if (c1->molkind != NULL)
                        if (c1->molkind[i1] != 0 && c1->MolID[i1] == c2->MolID[i2]) continue;  // TO-DO
                    */

                    unsigned atomkind2 = c2->atomkind[i2];

                    LJ_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, c1, i1, atomkind2, c2, i2,
                                                       PotEnergy, pottable);
                }
            }

            cell_t *c2 = c1;

            /* iterate over particles in cell c2=c1 with i2 < i1 */

            for (int i2=0; i2 < i1; i2++)
            {
                if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                /*
                if (c1->molkind != NULL)
                    if (c1->molkind[i1] != 0 && c1->MolID[i1] == c2->MolID[i2]) continue;  // TO-DO
                */

                unsigned atomkind2 = c2->atomkind[i2];

                LJ_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, c1, i1, atomkind2, c2, i2,
                                                   PotEnergy, pottable);
            }
        }
    }

    MPI_Allreduce(&PotEnergy, &md->GroupPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                            fmd_real_t sigma, fmd_real_t epsilon, fmd_real_t cutoff)
{
    LJ_6_12_t *lj = (LJ_6_12_t *)m_alloc(sizeof(LJ_6_12_t));
    lj->sig = sigma;
    lj->sig_sqr = sqrr(sigma);
    lj->eps = epsilon;
    lj->eps4 = 4. * epsilon;
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
