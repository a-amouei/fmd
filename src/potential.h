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
    POT_NONE,
    POT_LJ_6_12,
    POT_EAM_ALLOY,
    POT_MORSE
} potcat_t;

typedef struct
{
    potcat_t cat;
    void *data;
} fmd_pot_t;

typedef struct
{
    potcat_t cat;
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

typedef struct _list list_t;
typedef struct _bondkind bondkind_t;
typedef bondkind_t *bondkindp_t;
typedef struct _molkind molkind_t;

typedef struct
{
    unsigned atomkinds_num;
    atomkind_t *atomkinds;
    potpair_t **pottable;           // table of applied pots
    list_t *potcats;                // list of pot categories that are present in pottable
    unsigned potcats_num;
    list_t *potlist;                // list of all pots, whether applied or not
    fmd_bool_t hybridpasses[2];
    bondkindp_t *bondkinds;         // array of pointers to bondkinds that are defined
    unsigned bondkinds_num;         // size of 'bondkinds' array
    molkind_t *molkinds;
    unsigned  molkinds_num;
} potsys_t;

typedef struct _fmd fmd_t;

void fmd_potsys_free(fmd_t *md);
void fmd_potsys_init(fmd_t *md);
void fmd_pot_prepareForForceComp(fmd_t *md);
void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);

#endif /* POTENTIAL_H */
