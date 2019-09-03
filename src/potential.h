/*
  potential.h: This file is part of Free Molecular Dynamics

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

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "config.h"
#include "types.h"

typedef enum
{
    POTKIND_NONE,
    POTKIND_LJ_6_12,
    POTKIND_EAM_ALLOY,
    POTKIND_MORSE
} potkind_t;

typedef struct
{
    potkind_t kind;
    void *data;
} fmd_pot_t;

typedef struct
{
    potkind_t kind;
    void *data;
    unsigned iloc, jloc;    // local indexes inside the potential, used in potentials like EAM
} potpair_t;

typedef struct eam_element_t eam_element_t;

typedef struct
{
    double mass;
    fmd_string_t name;
    eam_element_t *eam_element;
} atomkind_t;

typedef struct list_t list_t;
typedef struct bond_t bond_t;

typedef struct
{
    unsigned atomkinds_num;
    atomkind_t *atomkinds;
    potpair_t **pottable;       // table of applied pots
    list_t *potkinds;           // list of pot kinds that are present in pottable
    unsigned potkinds_num;
    list_t *potlist;            // list of all pots, whether applied or not
    bond_t *bonds;              // array of bonds that are defined
    int hybridpasses[2];
} potsys_t;

typedef struct fmd_t fmd_t;

void fmd_potsys_free(fmd_t *md);
void fmd_potsys_init(fmd_t *md);
void fmd_pot_prepareForForceComp(fmd_t *md);
void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);

#endif /* POTENTIAL_H */
