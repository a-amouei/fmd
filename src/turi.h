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
    FMD_FIELD_MASS,
    FMD_FIELD_VCM,
    FMD_FIELD_TEMPERATURE,
    FMD_FIELD_NUMBER,
    FMD_FIELD_NUMBER_DENSITY,
    FMD_FIELD_TTM_TE,
    FMD_FIELD_TTM_XI
} fmd_field_t; /* category of the field */

typedef struct
{
    int field_index;
    fmd_field_t cat;
    int timestep;                    /* the last timestep this field was updated */
    fmd_real_t time;                 /* the last time this field was updated */
    unsigned dependcs_num;           /* number of the dependency fields */
    unsigned *dependcs;              /* indexes of the dependency fields */
    unsigned intervals_allreduce_num;
    fmd_real_t *intervals_allreduce; /* at these intervals a field update and "allreduce" MPI communication is done */
    fmd_bool_t allreduce_now;
    unsigned intervals_num;
    fmd_real_t *intervals;           /* intervals determine when to update the field and perform "reduce" MPI communication */
    fmd_array3D_t data;
    unsigned data_el_size;           /* size of each data element in bytes */
    datatype_t datatype;
} field_t;

typedef enum
{
    FMD_TURI_CUSTOM,
    FMD_TURI_TTM_TYPE1
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

typedef struct _ttm ttm_t;

typedef struct _turi
{
    int turi_index;
    fmd_turi_t cat;
    fmd_utuple_t tdims_global;  /* global dimenstions of the turi */
    fmd_utuple_t tdims_local;   /* dimenstions of the turi in current subdomain.
                                   turi-margin is included. */
    fmd_ituple_t tdims_local_nonmarg;  /* turi-margin is NOT included */


    fmd_ituple_t rank_of_lower_owner;  /* only if the current subdomain "owns" one or more turi-cells
                                          and cat == FMD_TURI_TTM, these three variables are initialized
                                          in a correct way. */
    fmd_ituple_t rank_of_upper_owner;
    fmd_btuple_t has_upper_lower_owner_procs;
    unsigned tcells_global_num;
    fmd_ituple_t itc_start;            /* width of turi-margin, refers to the first
                                          turi-cell in the interior of the local turi */
    fmd_ituple_t itc_stop;             /* first local index in the upper turi-margin */
    fmd_ituple_t itc_start_global;     /* global index of the first turi-cell in current subdomain */
    fmd_ituple_t itc_stop_global;      /* global index of the first turi-cell outside current subdomain */
    fmd_ituple_t itc_glob_to_loc;      /* itc_local[d] = itc_global[d] + itc_glob_to_loc[d] */
    fmd_ituple_t itc_start_owned;      /* local index of the first turi-cell "owned" by current subdomain */
    fmd_rtuple_t tcellh;               /* size of each cell of the turi (global) */
    fmd_real_t tcell_volume;           /* equals tcellh[0] x tcellh[1] x tcellh[2] */
    unsigned fields_num;
    field_t *fields;
    unsigned comms_num;         /* number of communicators = number of elements of the following array */
    turi_comm_t *comms;
    unsigned num_tcells_max;    /* number of turi-cells associated with the communicator having maximum number of turi-cells */
    turi_ownerscomm_t ownerscomm;
    fmd_real_t starttime;       /* do not update fields when time < starttime */
    fmd_real_t stoptime;        /* do not update fields when time > stoptime (unsless stoptime < starttime) */
    ttm_t *ttm;
} turi_t;

typedef struct _fmd fmd_t;

int _fmd_field_add(turi_t *t, fmd_field_t cat, fmd_real_t interval, fmd_bool_t allreduce);
void _fmd_turies_update(fmd_t *md, int time_iteration, fmd_real_t time,
                        fmd_bool_t Xupd, fmd_bool_t Vupd, fmd_bool_t Fupd);
void fmd_turi_free(fmd_t *md);

#endif /* TURI_H */
