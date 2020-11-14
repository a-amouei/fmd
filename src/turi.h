/*
  turi.h: This file is part of Free Molecular Dynamics

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

#ifndef TURI_H
#define TURI_H

#include "config.h"
#include "types.h"
#include "array.h"
#include <mpi.h>

typedef enum
{
    FIELD_DATA_UNSIGNED,
    FIELD_DATA_REAL,
    FIELD_DATA_RTUPLE
} field_data_type_t;

typedef enum
{
    FMD_FIELD_MASS,
    FMD_FIELD_VCM,
    FMD_FIELD_TEMPERATURE,
    FMD_FIELD_NUMBER
} fmd_field_t; /* category of the field */

typedef struct
{
    fmd_field_t cat;
    int timestep;               /* the last timestep this field was updated */
    unsigned dependcs_num;      /* number of the dependency fields */
    unsigned *dependcs;         /* indexes of the dependency fields */
    unsigned intervals_num;
    fmd_real_t *intervals;      /* intervals determine when to update the field */
    fmd_array3D_t data;
    unsigned data_el_size;      /* size of each data element in bytes */
    field_data_type_t data_type;
} field_t;

typedef enum
{
    FMD_TURI_CUSTOM,
    FMD_TURI_TTM
} fmd_turi_t; /* category of the turi */

typedef struct
{
    MPI_Comm comm;
    int commsize;
    int *pset;                  /* ranks of the processes of this comm in MD_comm */
    unsigned num_tcells;        /* number of turi-cells associated with this comm */
    fmd_ituple_t *itcs;         /* local indexes of turi-cells associated with this comm */
} turi_comm_t;

typedef struct
{
    fmd_ituple_t *owned_tcells;   /* list of local indexes of turi-cells "own"ed by current subdomain */
    int owned_tcells_num;         /* number of turi-cells "own"ed by current subdomain */
    MPI_Comm comm;
    int commsize;
    fmd_ituple_t *global_indexes; /* used when gattering data of all tcells (for root process) */
    int *recvcounts;              /* used when gattering data of all tcells (for root process) */
    int *displs;                  /* used when gattering data of all tcells (for root process) */
} turi_ownerscomm_t;

typedef struct _turi
{
    fmd_turi_t cat;
    fmd_ituple_t tdims_global;  /* global dimenstions of the turi */
    fmd_ituple_t tdims;         /* dimenstions of the turi in current subdomain.
                                   The subdomain may share some of its turi-cells with
                                   neighbor subdomains. */
    unsigned tcells_global_num;
    fmd_ituple_t tcell_start;   /* global index of the first turi-cell in current subdomain */
    fmd_ituple_t tcell_stop;    /* global index of the first turi-cell outside current subdomain */
    fmd_rtuple_t tcellh;        /* size of each cell of the turi (global) */
    unsigned fields_num;
    field_t *fields;
    unsigned comms_num;         /* number of communicators = number of elements of the following array */
    turi_comm_t *comms;
    unsigned num_tcells_max;    /* number of turi-cells associated with the communicator having maximum number of turi-cells */
    turi_ownerscomm_t ownerscomm;
    fmd_real_t starttime;       /* do not update fields when time < starttime */
    fmd_real_t stoptime;        /* do not update fields when time > stoptime (unsless stoptime < starttime) */
} turi_t;

typedef struct _fmd fmd_t;

fmd_handle_t fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz, fmd_real_t starttime, fmd_real_t stoptime);
fmd_handle_t fmd_field_add(fmd_t *md, fmd_handle_t turi, fmd_field_t cat, fmd_real_t interval);
void _fmd_turies_update(fmd_t *md);

#endif /* TURI_H */
