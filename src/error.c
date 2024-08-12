/*
  error.c: This file is part of Free Molecular Dynamics

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

#include <stdio.h>
#include "error.h"
#include "events.h"
#include "fmd-private.h"

#define FORMAT1 "FMD (%s, %s, %d): ERROR: "

static void error_set_common_members(
  fmd_event_params_error_t *err, bool major,
  fmd_string_t source, fmd_string_t func, int line)
{
    err->major = major;
    err->source = source;
    err->func = func;
    err->line = line;
    err->p1 = err->p2 = err->p3 = NULL;
}

void _fmd_error_unable_open_file(fmd_t *md, bool major, fmd_string_t source,
                                 fmd_string_t func, int line, fmd_string_t path)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Unable to open the file (%s)!\n", source, func, line, path);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_UNABLE_OPEN_FILE;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = path;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_unexpected_position(fmd_t *md, bool major, fmd_string_t source,
                                    fmd_string_t func, int line)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Unexpected particle position!\n", source, func, line);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_UNEXPECTED_POSITION;
    error_set_common_members(&err, major, source, func, line);

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_outside_real_interval(fmd_t *md, bool major, fmd_string_t source,
                                      fmd_string_t func, int line, fmd_string_t qname,
                                      fmd_real_t val, fmd_string_t intvname)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"The %s (%g) is outside the interval %s!\n",
                source, func, line, qname, val, intvname);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_OUTSIDE_REAL_INTERVAL;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = qname;
    err.p2 = &val;
    err.p3 = intvname;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_unacceptable_int_value(fmd_t *md, bool major, fmd_string_t source,
                                       fmd_string_t func, int line, fmd_string_t qname,
                                       int val)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Unacceptable %s (%d)!\n",
                source, func, line, qname, val);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_UNACCEPTABLE_INT_VALUE;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = qname;
    err.p2 = &val;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_file_corrupted(fmd_t *md, bool major, fmd_string_t source,
                               fmd_string_t func, int line, fmd_string_t ftype,
                               fmd_string_t path)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Not a healthy %s file (%s)!\n",
                source, func, line, ftype, path);
    
    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_FILE_CORRUPTED;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = ftype;
    err.p2 = path;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_unable_allocate_mem(fmd_t *md, bool major, fmd_string_t source,
                                    fmd_string_t func, int line, size_t size)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Unable to allocate memory (%zu bytes)!\n",
                source, func, line, size);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_UNABLE_ALLOCATE_MEM;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = &size;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_unsuccessful_hdf5(fmd_t *md, bool major, fmd_string_t source,
                                  fmd_string_t func, int line)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Unsuccessful HDF5 operation!\n", source, func, line);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_UNSUCCESSFUL_HDF5;
    error_set_common_members(&err, major, source, func, line);

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_function_failed(fmd_t *md, bool major, fmd_string_t source,
                                fmd_string_t func, int line, fmd_string_t name)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"The function %s() failed!\n", source, func, line, name);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_FUNCTION_FAILED;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = name;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_not_supported_yet(fmd_t *md, bool major, fmd_string_t source,
                                  fmd_string_t func, int line, fmd_string_t msg)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"%s\n", source, func, line, msg);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_NOT_SUPPORTED_YET;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = msg;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_unprepared(fmd_t *md, bool major, fmd_string_t source,
                           fmd_string_t func, int line, fmd_string_t name)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"Unprepared %s!\n", source, func, line, name);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_UNPREPARED;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = name;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}

void _fmd_error_wrong_potential(fmd_t *md, bool major, fmd_string_t source,
                                fmd_string_t func, int line, fmd_string_t pname,
                                int atomkind1, int atomkind2)
{
    if (md->ShowErrorMessages)
        fprintf(stderr, FORMAT1"The chosen %s potential does not support the (%d, %d) pair of atomkinds!\n",
                source, func, line, pname, atomkind1, atomkind2);

    if (md->EventHandler == NULL) return;

    fmd_event_params_error_t err;

    err.error = FMD_ERR_WRONG_POTENTIAL;
    error_set_common_members(&err, major, source, func, line);
    err.p1 = pname;
    err.p2 = &atomkind1;
    err.p3 = &atomkind2;

    md->EventHandler(md, FMD_EVENT_ERROR, md->userobject, (fmd_params_t *)&err);
}
