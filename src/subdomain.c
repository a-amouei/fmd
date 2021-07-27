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
#include "general.h"

fmd_real_t _fmd_convert_pos_to_subd_coord_1D(fmd_t *md, fmd_real_t pos, int d)
{
    fmd_real_t pos2, ref2, width1;

    width1 = (md->SubDomain.w[d] + 1) * md->cellh[d];
    ref2 = md->SubDomain.r[d] * width1;
    pos2 = pos - ref2;

    if (pos2 <= 0.0)
        return pos / width1;
    else
        return md->SubDomain.r[d] + pos2 / (md->SubDomain.w[d] * md->cellh[d]);
}

/* This function receives the position of a point in MD simulation box
   and determines its coordinates in subdomain space.
   For example if for a pos[] value we obtain s[] = {1.5, 0.5, 3.5}, it
   means that pos[] is placed right in the center of the subdomain with
   index of is[] = {1, 0, 3}. */
void _fmd_convert_pos_to_subd_coord(fmd_t *md, fmd_rtuple_t pos, fmd_rtuple_t s)
{
    for (int d=0; d<3; d++)
        s[d] = _fmd_convert_pos_to_subd_coord_1D(md, pos[d], d);
}

void fmd_subd_free(fmd_t *md)
{
    if (md->SubDomain.grid != NULL)
    {
        _fmd_cleanGridSegment(md->SubDomain.grid, _fmd_ThreeZeros_int, md->SubDomain.cell_num);
        _fmd_array_3d_free(&md->SubDomain.grid_array);
        md->SubDomain.grid = NULL;
    }
}

void fmd_subd_init(fmd_t *md)
{
    int d;

    // initialize is
    INDEX_3D(md->SubDomain.myrank, md->ns, md->SubDomain.is);
    // initialize rank_of_lower_subd and rank_of_upper_subd (neighbor processes)
    fmd_ituple_t istemp;
    for (d=0; d<3; d++)
        istemp[d] = md->SubDomain.is[d];
    for (d=0; d<3; d++)
    {
        istemp[d] = (md->SubDomain.is[d] - 1 + md->ns[d]) % md->ns[d];
        md->SubDomain.rank_of_lower_subd[d] = INDEX_FLAT(istemp, md->ns);
        istemp[d] = (md->SubDomain.is[d] + 1) % md->ns[d];
        md->SubDomain.rank_of_upper_subd[d] = INDEX_FLAT(istemp, md->ns);
        istemp[d] = md->SubDomain.is[d];
    }
    //
    for (d=0; d<3; d++)
    {
        int r, w;

        if (md->ns[d] == 1) md->SubDomain.ic_start[d] = 0; else md->SubDomain.ic_start[d] = 1;

        md->SubDomain.r[d] = r = md->nc[d] % md->ns[d];
        md->SubDomain.w[d] = w = md->nc[d] / md->ns[d];

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

    _fmd_array_3d_create(md->SubDomain.cell_num, sizeof(cell_t), DATATYPE_CELL, &md->SubDomain.grid_array);
    md->SubDomain.grid = (cell_t ***)md->SubDomain.grid_array.data;
    assert(md->SubDomain.grid != NULL);
    /* TO-DO: handle memory error */

    _fmd_initialize_grid(md->SubDomain.grid, md->SubDomain.cell_num[0], md->SubDomain.cell_num[1], md->SubDomain.cell_num[2]);
}
