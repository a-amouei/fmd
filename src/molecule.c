/*
  molecule.c: This file is part of Free Molecular Dynamics

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

#include "molecule.h"
#include "potential.h"
#include "base.h"

unsigned fmd_bond_addKind(fmd_t *md, fmd_bond_t cat, double coeffs[])
{
    unsigned i = md->potsys.bondkinds_num;

    md->potsys.bondkinds = (bondkindp_t *)realloc(md->potsys.bondkinds,
      (i+1) * sizeof(bondkindp_t));

    switch (cat)
    {
        case FMD_BOND_HARMONIC:
        {
            bondkind_harmonic_t *bkh = (bondkind_harmonic_t *)malloc(sizeof(bondkind_harmonic_t));
            bkh->cat = cat;
            bkh->k = coeffs[0];
            bkh->r0 = coeffs[1];
            md->potsys.bondkinds[i] = (bondkindp_t)bkh;
            break;
        }
    }

    md->potsys.bondkinds_num++;
    return i;
}

void fmd_bond_freeKinds(fmd_t *md)
{
    for (unsigned i=0; i < md->potsys.bondkinds_num; i++)
        free(md->potsys.bondkinds[i]);
    free(md->potsys.bondkinds);
    md->potsys.bondkinds_num = 0;
}

unsigned fmd_molecule_addKind(fmd_t *md, fmd_string_t name, unsigned AtomsNum, unsigned AtomKinds[])
{
    unsigned i = md->potsys.molkinds_num;

    md->potsys.molkinds = (molkind_t *)realloc(md->potsys.molkinds,
      (i+1) * sizeof(molkind_t));

    md->potsys.molkinds_num++;
    return i;
}
