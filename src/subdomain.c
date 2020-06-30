/*
  subdomain.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2020 Arham Amouye Foumani

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

#include "subdomain.h"
#include "base.h"

void fmd_subd_free(fmd_t *md)
{
    if (md->SubDomain.grid != NULL)
    {
        fmd_internal_freeGrid(md->SubDomain.grid, md->SubDomain.cell_num);
        md->SubDomain.grid = NULL;
    }
}

void fmd_subd_init(fmd_t *md)
{
    int d;

    // initialize is
    INVERSEINDEX(md->SubDomain.myrank, md->ns, md->SubDomain.is);
    // initialize rank_of_lower_subd and rank_of_upper_subd (neighbor processes)
    int istemp[3];
    for (d=0; d<3; d++)
        istemp[d] = md->SubDomain.is[d];
    for (d=0; d<3; d++)
    {
        istemp[d] = (md->SubDomain.is[d] - 1 + md->ns[d]) % md->ns[d];
        md->SubDomain.rank_of_lower_subd[d] = INDEX(istemp, md->ns);
        istemp[d] = (md->SubDomain.is[d] + 1) % md->ns[d];
        md->SubDomain.rank_of_upper_subd[d] = INDEX(istemp, md->ns);
        istemp[d] = md->SubDomain.is[d];
    }
    //
    for (d=0; d<3; d++)
    {
        int r, w;

        if (md->ns[d] == 1) md->SubDomain.ic_start[d] = 0; else md->SubDomain.ic_start[d] = 1;
        r = md->nc[d] % md->ns[d];
        w = md->nc[d] / md->ns[d];
        if (md->SubDomain.is[d] < r)
        {
            md->SubDomain.ic_stop[d] = md->SubDomain.ic_start[d] + w + 1;
            md->SubDomain.ic_global_firstcell[d] = md->SubDomain.is[d] * (w + 1);
        }
        else
        {
            md->SubDomain.ic_stop[d] = md->SubDomain.ic_start[d] + w;
            md->SubDomain.ic_global_firstcell[d] = md->SubDomain.is[d] * w + r;
        }
        md->SubDomain.cell_num[d] = md->SubDomain.ic_stop[d] + md->SubDomain.ic_start[d];
        md->SubDomain.cell_num_nonmarg[d] = md->SubDomain.ic_stop[d] - md->SubDomain.ic_start[d];
    }

    md->SubDomain.grid = fmd_internal_createGrid(md->SubDomain.cell_num);
}
