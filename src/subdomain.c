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

    width1 = (md->subd.w[d] + 1) * md->cellh[d];
    ref2 = md->subd.r[d] * width1;
    pos2 = pos - ref2;

    if (pos2 <= 0.0)
        return pos / width1;
    else
        return md->subd.r[d] + pos2 / (md->subd.w[d] * md->cellh[d]);
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

void _fmd_subd_free(fmd_t *md)
{
    if (md->subd.grid != NULL)
    {
        for (int ic=0; ic < md->subd.ncm; ic++)
            _fmd_cell_free(md->subd.grid + ic);

        free(md->subd.grid);
        md->subd.grid = NULL;

        _fmd_array_3d_free(&md->subd.gridp_array);
    }
}

static void init_cneighbs(Subdomain_t *s)
{
    fmd_ituple_t ic, jc;

    LOOP3D(ic, s->ic_start, s->ic_stop)
    {
        cell_t *c = ARRAY_ELEMENT(s->gridp, ic);
        int i=0;

        for (jc[0] = ic[0]-1; jc[0] <= ic[0]+1; jc[0]++)
            for (jc[1] = ic[1]-1; jc[1] <= ic[1]+1; jc[1]++)
                for (jc[2] = ic[2]-1; jc[2] <= ic[2]+1; jc[2]++)
                    c->cneighbs[i++] = ARRAY_ELEMENT(s->gridp, jc);
    }
}

static void init_gridp(Subdomain_t *s)
{
    fmd_ituple_t ic;
    unsigned imarg = s->nc, inonmarg = 0;

    LOOP3D(ic, _fmd_ThreeZeros_int, s->cell_num)
    {
        bool margin = false;

        for (int d=0; d<DIM; d++)
            if (ic[d] < s->ic_start[d] || ic[d] >= s->ic_stop[d])
            {
                margin = true;
                break;
            }

        if (margin)
            ARRAY_ELEMENT(s->gridp, ic) = s->grid + imarg++;
        else
        {
            ARRAY_ELEMENT(s->gridp, ic) = s->grid + inonmarg++;
        }
    }
}

void _fmd_subd_init(fmd_t *md)
{
    /* initialize is */
    INDEX_3D(md->subd.myrank, md->ns, md->subd.is);

    /* initialize rank_of_lower_subd and rank_of_upper_subd (neighbor processes) */
    fmd_ituple_t istemp;

    for (int d=0; d<DIM; d++)
        istemp[d] = md->subd.is[d];

    for (int d=0; d<DIM; d++)
    {
        if (!md->PBC[d] && (md->subd.is[d] == 0))
            md->subd.rank_of_lower_subd[d] = MPI_PROC_NULL;
        else
        {
            istemp[d] = (md->subd.is[d] - 1 + md->ns[d]) % md->ns[d];
            md->subd.rank_of_lower_subd[d] = INDEX_FLAT(istemp, md->ns);
        }

        if (!md->PBC[d] && (md->subd.is[d] == md->ns[d]-1))
            md->subd.rank_of_upper_subd[d] = MPI_PROC_NULL;
        else
        {
            istemp[d] = (md->subd.is[d] + 1) % md->ns[d];
            md->subd.rank_of_upper_subd[d] = INDEX_FLAT(istemp, md->ns);
        }

        istemp[d] = md->subd.is[d];
    }

    /* initialize other members of Subdomain_t */
    for (int d=0; d<DIM; d++)
    {
        int r, w;

        md->subd.ic_start[d] = 1;

        md->subd.r[d] = r = md->nc[d] % md->ns[d];
        md->subd.w[d] = w = md->nc[d] / md->ns[d];

        if (md->subd.is[d] < r)
        {
            md->subd.ic_stop[d] = md->subd.ic_start[d] + w + 1;
            md->subd.ic_global_firstcell[d] = md->subd.is[d] * (w + 1);
        }
        else
        {
            md->subd.ic_stop[d] = md->subd.ic_start[d] + w;
            md->subd.ic_global_firstcell[d] = md->subd.is[d] * w + r;
        }

        md->subd.cell_num[d] = md->subd.ic_stop[d] + md->subd.ic_start[d];
        md->subd.cell_num_nonmarg[d] = md->subd.ic_stop[d] - md->subd.ic_start[d];
    }

    md->subd.nc = md->subd.cell_num_nonmarg[0] *
                  md->subd.cell_num_nonmarg[1] *
                  md->subd.cell_num_nonmarg[2];

    md->subd.ncm = md->subd.cell_num[0] *
                   md->subd.cell_num[1] *
                   md->subd.cell_num[2];

    md->subd.NumberOfParticles = 0;

    /* create grid and gridp and initialize them */

    md->subd.grid = m_alloc(md->subd.ncm * sizeof(cell_t));

    for (int ic=0; ic < md->subd.ncm; ic++)
        _fmd_cell_init(&md->cellinfo, md->subd.grid + ic);

    _fmd_array_3d_create(md->subd.cell_num, sizeof(cell_t *), DATATYPE_CELLP, &md->subd.gridp_array);
    md->subd.gridp = (cell_t ****)md->subd.gridp_array.data;
    assert(md->subd.gridp != NULL); /* TO-DO: handle memory error */

    init_gridp(&md->subd);

    /* initialize cneighbs array of each cell */
    init_cneighbs(&md->subd);
}

void _fmd_conv_ic_loc_to_glob(fmd_t *md, fmd_ituple_t ic, fmd_ituple_t icglob)
{
    for (int d=0; d<DIM; d++)
        icglob[d] = ic[d] - md->subd.ic_start[d] + md->subd.ic_global_firstcell[d];
}
