/*
  cell.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2021 Arham Amouye Foumani

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

#include "cell.h"
#include "fmd-private.h"
#include "misc.h"
#include "types.h"
#include "general.h"

void _fmd_cell_create_force_arrays(cell_t *c, bool FembP_alter)
{
    if (c->F == NULL) c->F = m_alloc(c->capacity * sizeof(fmd_rtuple_t));

    if (FembP_alter)
    {
        if (c->vaream == NULL)
            c->vaream = m_alloc(c->capacity * sizeof(fmd_real_t));
        else
        {
            free(c->vaream);
            c->vaream = NULL;
        }
    }
}

static inline void realloc_arrays(cell_t *c)
{
    if (c->x != NULL) c->x = re_alloc(c->x, c->capacity * sizeof(fmd_rtuple_t));
    if (c->v != NULL) c->v = re_alloc(c->v, c->capacity * sizeof(fmd_rtuple_t));
    if (c->F != NULL) c->F = re_alloc(c->F, c->capacity * sizeof(fmd_rtuple_t));
    if (c->vaream != NULL) c->vaream = re_alloc(c->vaream, c->capacity * sizeof(fmd_real_t));
    if (c->GroupID != NULL) c->GroupID = re_alloc(c->GroupID, c->capacity * sizeof(int));
    if (c->AtomID != NULL) c->AtomID = re_alloc(c->AtomID, c->capacity * sizeof(unsigned));
    if (c->atomkind != NULL) c->atomkind = re_alloc(c->atomkind, c->capacity * sizeof(unsigned));
}

void _fmd_cell_init(cellinfo_t *cinfo, cell_t *c)
{
    c->x = (cinfo->x_active ? m_alloc(sizeof(fmd_rtuple_t)) : NULL);
    c->v = (cinfo->v_active ? m_alloc(sizeof(fmd_rtuple_t)) : NULL);
    c->F = (cinfo->F_active ? m_alloc(sizeof(fmd_rtuple_t)) : NULL);
    c->vaream = (cinfo->vaream_active ? m_alloc(sizeof(fmd_real_t)) : NULL);
    c->GroupID = (cinfo->GroupID_active ? m_alloc(sizeof(int)) : NULL);
    c->AtomID = (cinfo->AtomID_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->atomkind = (cinfo->atomkind_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->capacity = 1;
    c->parts_num = 0;
}

void _fmd_cellinfo_init(cellinfo_t *cinfo)
{
    cinfo->x_active = true;
    cinfo->v_active = true;
    cinfo->F_active = false;
    cinfo->vaream_active = false;
    cinfo->GroupID_active = true;
    cinfo->AtomID_active = true;
    cinfo->atomkind_active = true;
}

void _fmd_cell_resize(fmd_t *md, cell_t *c)
{
    c->capacity = c->parts_num + md->cell_increment;

    realloc_arrays(c);
}

void _fmd_cell_resize_exact(cell_t *c)
{
    c->capacity = c->parts_num;

    realloc_arrays(c);
}

void _fmd_cell_minimize(fmd_t *md, cell_t *c)
{
    c->parts_num = 0;

    if (c->capacity > md->cell_increment)
    {
        c->capacity = md->cell_increment;
        realloc_arrays(c);
    }
}

void _fmd_cell_free(cell_t *c)
{
    c->parts_num = 0;
    c->capacity = 0;

    free(c->x);
    c->x = NULL;

    free(c->v);
    c->v = NULL;

    free(c->F);
    c->F = NULL;

    free(c->vaream);
    c->vaream = NULL;

    free(c->GroupID);
    c->GroupID = NULL;

    free(c->AtomID);
    c->AtomID = NULL;

    free(c->atomkind);
    c->atomkind = NULL;
}

void _fmd_cell_copy_atom_from_cell_to_cell(cell_t *cfrom, unsigned ifrom, cell_t *cto, unsigned ito)
{
    if (cfrom->x != NULL)
        for (int d=0; d<DIM; d++)
            POS(cto, ito, d) = POS(cfrom, ifrom, d);

    if (cfrom->v != NULL)
        for (int d=0; d<DIM; d++)
            VEL(cto, ito, d) = VEL(cfrom, ifrom, d);

    if (cfrom->F != NULL)
        for (int d=0; d<DIM; d++)
            FRC(cto, ito, d) = FRC(cfrom, ifrom, d);

    if (cfrom->vaream != NULL) cto->vaream[ito] = cfrom->vaream[ifrom];
    if (cfrom->GroupID != NULL) cto->GroupID[ito] = cfrom->GroupID[ifrom];
    if (cfrom->AtomID != NULL) cto->AtomID[ito] = cfrom->AtomID[ifrom];
    if (cfrom->atomkind != NULL) cto->atomkind[ito] = cfrom->atomkind[ifrom];
}

void _fmd_cell_remove_atom(fmd_t *md, cell_t *c, unsigned ind)
{
    c->parts_num--;

    if (c->parts_num != ind)
        _fmd_cell_copy_atom_from_cell_to_cell(c, c->parts_num, c, ind);

    if (c->parts_num + md->cell_increment < c->capacity)
        _fmd_cell_resize(md, c);
}
