/*
  fmd-private.h: This file is part of Free Molecular Dynamics

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

#ifndef FMD_PRIVATE_H
#define FMD_PRIVATE_H

#include "config.h"

#include <mpi.h>
#include <stdio.h>
#include "types.h"
#include "subdomain.h"
#include "potential.h"
#include "cell.h"
#include "events.h"
#include "h5.h"

#define MAX_PATH_LENGTH                     256
#define FMD_GROUP_ALL                        -1

typedef enum
{
    FMD_SCM_XYZ_PARTICLESNUM,
    FMD_SCM_XYZ_SEPARATE,
    FMD_SCM_CSV,
    FMD_SCM_VTF
} fmd_SaveConfigMode_t;

typedef struct _fmd_timer fmd_timer_t;

typedef struct
{
    MPI_Datatype mpi_ituple;
    MPI_Datatype mpi_rtuple;
    MPI_Datatype mpi_configa;
    MPI_Datatype mpi_statea;
} mpi_types_t;

typedef struct _turi turi_t;

struct _fmd
{
    Subdomain_t Subdomain;
    potsys_t potsys;
    cell_t ***global_grid;
    cellinfo_t cellinfo;
    fmd_EventHandler_t EventHandler;
    mpi_types_t mpi_types;
    unsigned timers_num;
    fmd_timer_t *timers;
    unsigned turies_num;
    turi_t *turies;
    turi_t *active_ttm_turi;
    fmd_bool_t GlobalGridExists;
    fmd_bool_t BoxSizeDetermined;
    fmd_bool_t PBCdetermined;
    fmd_real_t CutoffRadius;
    fmd_real_t time;
    int time_iteration;
    fmd_real_t timestep;
    unsigned TotalNoOfParticles;
    unsigned TotalNoOfMolecules;
    fmd_real_t GroupTemperature;          /* for the "active group" only */
    fmd_bool_t Is_MD_process;
    fmd_bool_t Is_MD_comm_root;
    int LOP_iteration;                    /* must be initialized with zero */
    int LOP_period;
    MPI_Comm MD_comm;
    fmd_real_t GroupKineticEnergy;
    fmd_real_t GroupPotentialEnergy;
    int world_rank;
    int world_numprocs;
    fmd_real_t DesiredTemperature;
    fmd_ituple_t PBC;
    fmd_ituple_t ns;                      // number of subdomains = ns[0] x ns[1] x ns[2]
    fmd_rtuple_t l;                       // size of the simulation box
    fmd_ituple_t nc;                      // number of grid cells in the simulation box
    fmd_rtuple_t cellh;                   // size of one single grid cell
    fmd_bool_t UseAutoStep;
    fmd_real_t AutoStepSensitivity;
    char SaveDirectory[MAX_PATH_LENGTH];
    fmd_real_t BerendsenThermostatParam;
    fmd_SaveConfigMode_t SaveConfigMode;
    FILE *ConfigFilep;
    fmd_real_t WallTimeOrigin;
    int ActiveGroup;
    unsigned GroupParticlesNum;
    fmd_rtuple_t GroupMomentum;
    fmd_bool_t ParticlesDistributed;
    fmd_bool_t MPI_initialized_by_me;
    h5_dataspaces_t h5_dataspaces;
    int cell_increment;
    int _OldNumberOfParticles;
    int _FileIndex;
};

typedef struct _fmd fmd_t;

#endif /* FMD_PRIVATE_H */
