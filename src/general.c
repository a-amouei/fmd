/*
  general.c: This file is part of Free Molecular Dynamics

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

#include <mpi.h>
#include "general.h"

const fmd_itriple_t _fmd_ThreeZeros_int = {0, 0, 0};

const fmd_ftriple_t _fmd_ThreeZeros_float = {0.0, 0.0, 0.0};

FILE *f_open(char *filename, char *modes)
{
    FILE *fp = fopen(filename, modes);

    if (fp == NULL)
    {
        fprintf(stderr, "ERROR: Unable to open %s!\n", filename);
        MPI_Abort(MPI_COMM_WORLD, ERROR_UNABLE_OPEN_FILE);
    }

    return fp;
}
