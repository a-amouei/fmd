/*
  misc.h: This file is part of Free Molecular Dynamics

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

#ifndef MISC_H
#define MISC_H

#include "config.h"

/* Global macroes and symbolic constants */

#define LIGHT_SPEED                 2.9979245800e+06   // (ang / ps)
#define METER_PER_SECOND            1e-2               // (ang / ps)
#define JOULE_PER_METER2            6.2415091259e-02   // (eV / ang^2)
#define JOULE_PER_METER3_KELVIN     6.2415091259e-12   // (eV / ang^3 Kelvin)
#define PASCAL                      6.2415091259e-12   // (mass_unit / ang ps^2)
#define MD_CHARGE_UNIT              1.2657711566e-10   // x (electrostatic unit of charge (esu) = statcoulomb)
#define E_CHARGE                    3.7946864629e+00   // electron charge in MD electric charge unit
#define E_MASS                      5.6856300621e-08   // electron mass in MD mass unit
#define HBAR                        6.5821195136e-04   // Planck constant devided by 2*pi in MD units (eV ps)
#define GRAM_PER_CM3                6.2415091259e-05   // x (mass_unit / ang^3)
#define NEWTON_PER_METER            6.2415091259e-02   // x (eV / ang^2)

/* */

typedef struct _cell cell_t;
typedef struct _cellinfo cellinfo_t;
typedef struct _fmd fmd_t;

/* Functions */

void _fmd_createGlobalGrid(fmd_t *md);
void _fmd_refreshGrid(fmd_t *md);

#endif /* MISC_H */
