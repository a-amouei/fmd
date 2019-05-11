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

typedef char fmd_atomkind_name_t[17];

typedef struct
{
    double eps;
    double sig;
    double cutoff_sqr;
} LJ_6_12_t;

typedef struct
{
    double mass;
    fmd_atomkind_name_t name;
} atomkind_t;

typedef struct
{
    unsigned atomkinds_num;
    atomkind_t *atomkinds;
    LJ_6_12_t **lj_6_12;
} potential_t;

typedef struct fmd_sys_t fmd_sys_t;

void fmd_pot_free(fmd_sys_t *sysp);
void fmd_pot_init(fmd_sys_t *sysp);

#endif /* POTENTIAL_H */
