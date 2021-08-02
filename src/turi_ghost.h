/*
  turi_ghost.h: This file is part of Free Molecular Dynamics

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

#ifndef TURI_GHOST_H
#define TURI_GHOST_H

#include <stddef.h>
#include "config.h"
#include "types.h"

typedef struct _turi turi_t;
typedef struct _fmd fmd_t;

typedef void (*fields_packer_t)(fmd_t *md, turi_t *t, fmd_ituple_t vitc_start, fmd_ituple_t vitc_stop,
                                size_t *size, fmd_pointer_t *out);
typedef void (*fields_unpacker_t)(fmd_t *md, turi_t *t, fmd_ituple_t vitc_start, fmd_ituple_t vitc_stop, fmd_pointer_t in);

void _fmd_turi_update_ghosts(fmd_t *md, turi_t *t, fields_packer_t pack, fields_unpacker_t unpack);

#endif /* TURI_GHOST_H */
