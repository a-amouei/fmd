/*
  cell.h: This file is part of Free Molecular Dynamics

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

#ifndef CELL_H
#define CELL_H

#include "config.h"
#include "types.h"

#define COMP(i, d) ((i)*DIM + (d))
#define FRC(cell,i,d) ((cell)->F[COMP((i), (d))])
#define POS(cell,i,d) ((cell)->x[COMP((i), (d))])
#define VEL(cell,i,d) ((cell)->v[COMP((i), (d))])

typedef struct _list list_t;

struct _cell
{
    unsigned parts_num;
    unsigned capacity;
    fmd_real_t *x, *v, *F;
    fmd_real_t *FembPrime;
    int *GroupID;
    unsigned *AtomID;
    unsigned *atomkind;
    unsigned *molkind;
    unsigned *MolID;
    unsigned *AtomIDlocal;    /* Atom ID in molecule */
    list_t **neighbors;       /* for each molecule atom in cell there is a list of neighbors;
                                 each data pointer in this list points to a mol_atom_neighbor_t */
};

typedef struct _cell cell_t;

struct _cellinfo
{
    bool x_active;
    bool v_active;
    bool F_active;
    bool FembPrime_active;
    bool GroupID_active;
    bool AtomID_active;
    bool atomkind_active;
    bool molkind_active;
    bool MolID_active;
    bool AtomIDlocal_active;
    bool neighbors_active;
};

typedef struct _cellinfo cellinfo_t;

typedef struct _fmd fmd_t;

void _fmd_cell_resize(fmd_t *md, cell_t *c);
void _fmd_cell_resize_exact(cell_t *c);
void _fmd_cell_minimize(fmd_t *md, cell_t *c);
void _fmd_cell_free(cell_t *c);
void _fmd_cell_init(cellinfo_t *cinfo, cell_t *c);
void _fmd_cellinfo_init(cellinfo_t *cinfo);
void _fmd_cell_remove_atom(fmd_t *md, cell_t *c, unsigned ind);
void _fmd_cell_copy_atom_from_cell_to_cell(cell_t *cfrom, unsigned ifrom, cell_t *cto, unsigned ito);
void _fmd_cell_create_force_arrays(cell_t *c, bool FembP_alter);

inline unsigned _fmd_cell_new_particle(fmd_t *md, cell_t *c)
{
    if (c->parts_num == c->capacity) _fmd_cell_resize(md, c);

    return c->parts_num++;
}

#endif /* CELL_H */
