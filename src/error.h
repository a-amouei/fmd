/*
  error.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2024 Arham Amouye Foumani

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

#ifndef ERROR_H
#define ERROR_H

#include "config.h"
#include "types.h"

typedef enum
{
    FMD_ERR_UNEXPECTED_POSITION = 1,
    FMD_ERR_UNABLE_OPEN_FILE,
    FMD_ERR_UNABLE_ALLOCATE_MEM,
    FMD_ERR_FILE_CORRUPTED,
    FMD_ERR_OUTSIDE_REAL_INTERVAL,
    FMD_ERR_UNACCEPTABLE_INT_VALUE,
    FMD_ERR_UNSUCCESSFUL_HDF5,
    FMD_ERR_FUNCTION_FAILED,
    FMD_ERR_NOT_SUPPORTED_YET,
    FMD_ERR_UNPREPARED,
    FMD_ERR_WRONG_POTENTIAL
} fmd_error_t;

typedef struct _fmd fmd_t;

void _fmd_error_unable_open_file(fmd_t *md, bool major, fmd_string_t source,
                                 fmd_string_t func, int line, fmd_string_t path);

void _fmd_error_unexpected_position(fmd_t *md, bool major, fmd_string_t source,
                                    fmd_string_t func, int line);

void _fmd_error_outside_real_interval(fmd_t *md, bool major, fmd_string_t source,
                                      fmd_string_t func, int line, fmd_string_t qname,
                                      fmd_real_t val, fmd_string_t intvname);

void _fmd_error_unacceptable_int_value(fmd_t *md, bool major, fmd_string_t source,
                                       fmd_string_t func, int line, fmd_string_t qname,
                                       int val);

void _fmd_error_file_corrupted(fmd_t *md, bool major, fmd_string_t source,
                               fmd_string_t func, int line, fmd_string_t ftype,
                               fmd_string_t path);

void _fmd_error_unable_allocate_mem(fmd_t *md, bool major, fmd_string_t source,
                                    fmd_string_t func, int line, size_t size);

void _fmd_error_unsuccessful_hdf5(fmd_t *md, bool major, fmd_string_t source,
                                  fmd_string_t func, int line);

void _fmd_error_function_failed(fmd_t *md, bool major, fmd_string_t source,
                                fmd_string_t func, int line, fmd_string_t name);

void _fmd_error_not_supported_yet(fmd_t *md, bool major, fmd_string_t source,
                                  fmd_string_t func, int line, fmd_string_t msg);

void _fmd_error_unprepared(fmd_t *md, bool major, fmd_string_t source,
                           fmd_string_t func, int line, fmd_string_t name);

void _fmd_error_wrong_potential(fmd_t *md, bool major, fmd_string_t source,
                                fmd_string_t func, int line, fmd_string_t pname,
                                int atomkind1, int atomkind2);

#endif /* ERROR_H */
