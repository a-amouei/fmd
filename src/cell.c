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
#include "base.h"
#include "types.h"
#include "general.h"

static inline void realloc_arrays(cell_t *c)
{
    if (c->x != NULL) c->x = re_alloc(c->x, c->capacity * sizeof(fmd_rtuple_t));
    if (c->v != NULL) c->v = re_alloc(c->v, c->capacity * sizeof(fmd_rtuple_t));
    if (c->F != NULL) c->F = re_alloc(c->F, c->capacity * sizeof(fmd_rtuple_t));
    if (c->FembPrime != NULL) c->FembPrime = re_alloc(c->FembPrime, c->capacity * sizeof(fmd_real_t));
    if (c->GroupID != NULL) c->GroupID = re_alloc(c->GroupID, c->capacity * sizeof(int));
    if (c->AtomID != NULL) c->AtomID = re_alloc(c->AtomID, c->capacity * sizeof(unsigned));
    if (c->atomkind != NULL) c->atomkind = re_alloc(c->atomkind, c->capacity * sizeof(unsigned));
    if (c->molkind != NULL)
    {
        c->molkind = re_alloc(c->molkind, c->capacity * sizeof(unsigned));
        c->MolID = re_alloc(c->MolID, c->capacity * sizeof(unsigned));
        c->AtomIDlocal = re_alloc(c->AtomIDlocal, c->capacity * sizeof(unsigned));
        c->neighbors = re_alloc(c->neighbors, c->capacity * sizeof(list_t *));
    }
}

void _fmd_cell_init(cellinfo_t *cinfo, cell_t *c)
{
    c->x = (cinfo->x_active ? m_alloc(sizeof(fmd_rtuple_t)) : NULL);
    c->v = (cinfo->v_active ? m_alloc(sizeof(fmd_rtuple_t)) : NULL);
    c->F = (cinfo->F_active ? m_alloc(sizeof(fmd_rtuple_t)) : NULL);
    c->FembPrime = (cinfo->FembPrime_active ? m_alloc(sizeof(fmd_real_t)) : NULL);
    c->GroupID = (cinfo->GroupID_active ? m_alloc(sizeof(int)) : NULL);
    c->AtomID = (cinfo->AtomID_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->atomkind = (cinfo->atomkind_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->molkind = (cinfo->molkind_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->MolID = (cinfo->MolID_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->AtomIDlocal = (cinfo->AtomIDlocal_active ? m_alloc(sizeof(unsigned)) : NULL);
    c->neighbors = (cinfo->neighbors_active ? m_alloc(sizeof(list_t *)) : NULL);
    c->capacity = 1;
    c->parts_num = 0;
}

void _fmd_cellinfo_init(cellinfo_t *cinfo)
{
    cinfo->x_active = FMD_TRUE;
    cinfo->v_active = FMD_TRUE;
    cinfo->F_active = FMD_TRUE;
    cinfo->FembPrime_active = FMD_TRUE;
    cinfo->GroupID_active = FMD_TRUE;
    cinfo->AtomID_active = FMD_TRUE;
    cinfo->atomkind_active = FMD_TRUE;
    cinfo->molkind_active = FMD_FALSE;
    cinfo->MolID_active = FMD_FALSE;
    cinfo->AtomIDlocal_active = FMD_FALSE;
    cinfo->neighbors_active = FMD_FALSE;
}

void _fmd_cell_resize(fmd_t *md, cell_t *c)
{
    c->capacity = c->parts_num + md->cell_increment;

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

    free(c->FembPrime);
    c->FembPrime = NULL;

    free(c->GroupID);
    c->GroupID = NULL;

    free(c->AtomID);
    c->AtomID = NULL;

    free(c->atomkind);
    c->atomkind = NULL;

    free(c->molkind);
    c->molkind = NULL;

    free(c->MolID);
    c->MolID = NULL;

    free(c->AtomIDlocal);
    c->AtomIDlocal = NULL;

    free(c->neighbors);
    c->neighbors = NULL;
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

    if (cfrom->FembPrime != NULL) cto->FembPrime[ito] = cfrom->FembPrime[ifrom];
    if (cfrom->GroupID != NULL) cto->GroupID[ito] = cfrom->GroupID[ifrom];
    if (cfrom->AtomID != NULL) cto->AtomID[ito] = cfrom->AtomID[ifrom];
    if (cfrom->atomkind != NULL) cto->atomkind[ito] = cfrom->atomkind[ifrom];

    if (cfrom->molkind != NULL)
    {
        cto->molkind[ito] = cfrom->molkind[ifrom];
        cto->MolID[ito] = cfrom->MolID[ifrom];
        cto->AtomIDlocal[ito] = cfrom->AtomIDlocal[ifrom];
        cto->neighbors[ito] = cfrom->neighbors[ifrom];
    }
}

void _fmd_cell_remove_atom(fmd_t *md, cell_t *c, unsigned ind)
{
    c->parts_num--;

    if (c->parts_num != ind)
        _fmd_cell_copy_atom_from_cell_to_cell(c, c->parts_num, c, ind);

    if (c->parts_num + md->cell_increment < c->capacity)
        _fmd_cell_resize(md, c);
}
