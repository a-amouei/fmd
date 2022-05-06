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

static void clean_grid_cells(cell_t ***grid, fmd_utuple_t ex)
{
    fmd_ituple_t ic;

    LOOP3D(ic, _fmd_ThreeZeros_int, ex)
        _fmd_cell_free(&grid[ic[0]][ic[1]][ic[2]]);
}

void fmd_subd_free(fmd_t *md)
{
    if (md->SubDomain.grid != NULL)
    {
        clean_grid_cells(md->SubDomain.grid, md->SubDomain.cell_num);
        _fmd_array_3d_free(&md->SubDomain.grid_array);
        md->SubDomain.grid = NULL;
    }
}

void fmd_subd_init(fmd_t *md)
{
    /* initialize is */
    INDEX_3D(md->SubDomain.myrank, md->ns, md->SubDomain.is);

    /* initialize rank_of_lower_subd and rank_of_upper_subd (neighbor processes) */
    fmd_ituple_t istemp;

    for (int d=0; d<DIM; d++)
        istemp[d] = md->SubDomain.is[d];

    for (int d=0; d<DIM; d++)
    {
        if (!md->PBC[d] && (md->SubDomain.is[d] == 0))
            md->SubDomain.rank_of_lower_subd[d] = MPI_PROC_NULL;
        else
        {
            istemp[d] = (md->SubDomain.is[d] - 1 + md->ns[d]) % md->ns[d];
            md->SubDomain.rank_of_lower_subd[d] = INDEX_FLAT(istemp, md->ns);
        }

        if (!md->PBC[d] && (md->SubDomain.is[d] == md->ns[d]-1))
            md->SubDomain.rank_of_upper_subd[d] = MPI_PROC_NULL;
        else
        {
            istemp[d] = (md->SubDomain.is[d] + 1) % md->ns[d];
            md->SubDomain.rank_of_upper_subd[d] = INDEX_FLAT(istemp, md->ns);
        }

        istemp[d] = md->SubDomain.is[d];
    }

    /*  */
    for (int d=0; d<DIM; d++)
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

    _fmd_initialize_grid(md->SubDomain.grid, &md->cellinfo, md->SubDomain.cell_num[0],
                         md->SubDomain.cell_num[1], md->SubDomain.cell_num[2]);

    md->SubDomain.NumberOfParticles = 0;
}

void _fmd_conv_ic_loc_to_glob(fmd_t *md, fmd_ituple_t ic, fmd_ituple_t icglob)
{
    for (int d=0; d<3; d++)
        icglob[d] = ic[d] - md->SubDomain.ic_start[d] + md->SubDomain.ic_global_firstcell[d];
}
