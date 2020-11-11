/*
  molecule.h: This file is part of Free Molecular Dynamics

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

#ifndef MOLECULE_H
#define MOLECULE_H

#include "config.h"
#include "types.h"

typedef enum
{
    FMD_BOND_HARMONIC
} fmd_bond_t;   // category of a bond

typedef struct _bondkind
{
    fmd_bond_t cat;
} bondkind_t;

typedef bondkind_t *bondkindp_t;

typedef struct
{
    fmd_bond_t cat;
    fmd_real_t k;
    fmd_real_t r0;
} bondkind_harmonic_t;

typedef struct _list list_t;

typedef struct
{
    unsigned LocalID;
    unsigned atomkind;
    fmd_rtuple_t position;
    list_t *neighbors;   // each data pointer in this list points to a molkind_atom_neighbor_t
} molkind_atom_t;

typedef struct
{
    molkind_atom_t *atom;
    bondkind_t *bond;
} molkind_atom_neighbor_t;

typedef struct _molkind
{
    unsigned atoms_num;
    fmd_string_t name;
    unsigned **distances;
    molkind_atom_t *atoms;
} molkind_t;

typedef struct _fmd fmd_t;

void fmd_bond_freeKinds(fmd_t *md);
void _fmd_matt_updateAtomNeighbors(fmd_t *md);
void fmd_dync_computeBondForce(fmd_t *md);

#endif /* MOLECULE_H */
