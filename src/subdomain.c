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
#include "fmd-private.h"
#include "misc.h"
#include "general.h"

fmd_real_t _fmd_convert_pos_to_subd_coord_1D(fmd_t *md, fmd_real_t pos, int d)
{
    fmd_real_t pos2, ref2, width1;

    width1 = (md->Subdomain.w[d] + 1) * md->cellh[d];
    ref2 = md->Subdomain.r[d] * width1;
    pos2 = pos - ref2;

    if (pos2 <= 0.0)
        return pos / width1;
    else
        return md->Subdomain.r[d] + pos2 / (md->Subdomain.w[d] * md->cellh[d]);
}

/* This function receives the position of a point in MD simulation box
   and determines its coordinates in subdomain space.
   For example if for a pos[] value we obtain s[] = {1.5, 0.5, 3.5}, it
   means that pos[] is placed right in the center of the subdomain with
   index of is[] = {1, 0, 3}. */
void _fmd_convert_pos_to_subd_coord(fmd_t *md, fmd_rtuple_t pos, fmd_rtuple_t s)
{
    for (int d=0; d<DIM; d++)
        s[d] = _fmd_convert_pos_to_subd_coord_1D(md, pos[d], d);
}

static void clean_grid_cells(cell_t ***grid, fmd_utuple_t ex)
{
    fmd_ituple_t ic;

    LOOP3D(ic, _fmd_ThreeZeros_int, ex)
        _fmd_cell_free(&ARRAY_ELEMENT(grid, ic));
}

void _fmd_subd_free(fmd_t *md)
{
    if (md->Subdomain.grid != NULL)
    {
        clean_grid_cells(md->Subdomain.grid, md->Subdomain.cell_num);
        _fmd_array_3d_free(&md->Subdomain.grid_array);
        md->Subdomain.grid = NULL;
    }
}

void _fmd_subd_init(fmd_t *md)
{
    /* initialize is */
    INDEX_3D(md->Subdomain.myrank, md->ns, md->Subdomain.is);

    /* initialize rank_of_lower_subd and rank_of_upper_subd (neighbor processes) */
    fmd_ituple_t istemp;

    for (int d=0; d<DIM; d++)
        istemp[d] = md->Subdomain.is[d];

    for (int d=0; d<DIM; d++)
    {
        if (!md->PBC[d] && (md->Subdomain.is[d] == 0))
            md->Subdomain.rank_of_lower_subd[d] = MPI_PROC_NULL;
        else
        {
            istemp[d] = (md->Subdomain.is[d] - 1 + md->ns[d]) % md->ns[d];
            md->Subdomain.rank_of_lower_subd[d] = INDEX_FLAT(istemp, md->ns);
        }

        if (!md->PBC[d] && (md->Subdomain.is[d] == md->ns[d]-1))
            md->Subdomain.rank_of_upper_subd[d] = MPI_PROC_NULL;
        else
        {
            istemp[d] = (md->Subdomain.is[d] + 1) % md->ns[d];
            md->Subdomain.rank_of_upper_subd[d] = INDEX_FLAT(istemp, md->ns);
        }

        istemp[d] = md->Subdomain.is[d];
    }

    /*  */
    for (int d=0; d<DIM; d++)
    {
        int r, w;

        if (md->ns[d] == 1) md->Subdomain.ic_start[d] = 0; else md->Subdomain.ic_start[d] = 1;

        md->Subdomain.r[d] = r = md->nc[d] % md->ns[d];
        md->Subdomain.w[d] = w = md->nc[d] / md->ns[d];

        if (md->Subdomain.is[d] < r)
        {
            md->Subdomain.ic_stop[d] = md->Subdomain.ic_start[d] + w + 1;
            md->Subdomain.ic_global_firstcell[d] = md->Subdomain.is[d] * (w + 1);
        }
        else
        {
            md->Subdomain.ic_stop[d] = md->Subdomain.ic_start[d] + w;
            md->Subdomain.ic_global_firstcell[d] = md->Subdomain.is[d] * w + r;
        }

        md->Subdomain.cell_num[d] = md->Subdomain.ic_stop[d] + md->Subdomain.ic_start[d];
        md->Subdomain.cell_num_nonmarg[d] = md->Subdomain.ic_stop[d] - md->Subdomain.ic_start[d];
    }

    _fmd_array_3d_create(md->Subdomain.cell_num, sizeof(cell_t), DATATYPE_CELL, &md->Subdomain.grid_array);
    md->Subdomain.grid = (cell_t ***)md->Subdomain.grid_array.data;
    assert(md->Subdomain.grid != NULL);
    /* TO-DO: handle memory error */

    _fmd_initialize_grid(md->Subdomain.grid, &md->cellinfo, md->Subdomain.cell_num[0],
                         md->Subdomain.cell_num[1], md->Subdomain.cell_num[2]);

    md->Subdomain.NumberOfParticles = 0;
}

void _fmd_conv_ic_loc_to_glob(fmd_t *md, fmd_ituple_t ic, fmd_ituple_t icglob)
{
    for (int d=0; d<3; d++)
        icglob[d] = ic[d] - md->Subdomain.ic_start[d] + md->Subdomain.ic_global_firstcell[d];
}
