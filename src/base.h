/*
  base.h: This file is part of Free Molecular Dynamics

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

#ifndef BASE_H
#define BASE_H

#include "config.h"

#include <tgmath.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <limits.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "potential.h"
#include "subdomain.h"
#include "events.h"
#include "h5.h"
#include "cell.h"

/* Global macroes and symbolic constants */

#define K_BOLTZMANN                 8.6173303e-5       // (eV / Kelvin)
#define LIGHT_SPEED                 2.9979245800e+06   // (ang / ps)
#define METER_PER_SECOND            1e-2               // (ang / ps)
#define JOULE_PER_METER2            6.2415091259e-02   // (eV / ang^2)
#define JOULE_PER_METER3_KELVIN     6.2415091259e-12   // (eV / ang^3 Kelvin)
#define MD_MASS_UNIT                9.6485332907e+03   // x unified atomic mass unit
#define PASCAL                      6.2415091259e-12   // (mass_unit / ang ps^2)
#define MD_CHARGE_UNIT              1.2657711566e-10   // x (electrostatic unit of charge (esu) = statcoulomb)
#define E_CHARGE                    3.7946864629e+00   // electron charge in MD electric charge unit
#define E_MASS                      5.6856300621e-08   // electron mass in MD mass unit
#define HBAR                        6.5821195136e-04   // Planck constant devided by 2*pi in MD units (eV ps)
#define GRAM_PER_CM3                6.2415091259e-05   // x (mass_unit / ang^3)
#define NEWTON_PER_METER            6.2415091259e-02   // x (eV / ang^2)
#define MAX_PATH_LENGTH             256

#define ACTIVE_GROUP_ALL                        -1

/* typedefs & structs */

/*typedef struct
{
    //fmd_rtuple_t x_bak;
    //fmd_rtuple_t v_bak;
    float LocOrdParam;
    float x_avgd[3];
    unsigned AtomID;
    unsigned atomkind;
    unsigned molkind;
    unsigned MolID;
    unsigned AtomID_local;
} particle_core_t;

typedef struct _particle
{
    particle_core_t core;
    fmd_real_t FembPrime;
    float LocOrdParamAvg;
    list_t *neighbors;       *//* each data pointer in this list points to a mol_atom_neighbor_t */
/*} particle_t;*/

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
    SubDomain_t SubDomain;
    potsys_t potsys;
    cell_t ***global_grid;
    cellinfo_t cellinfo;
    fmd_EventHandler_t EventHandler;
    mpi_types_t mpi_types;
    unsigned timers_num;
    fmd_timer_t *timers;
    unsigned turies_num;
    turi_t *turies;
    fmd_bool_t GlobalGridExists;
    fmd_bool_t BoxSizeDetermined;
    fmd_bool_t PBCdetermined;
    fmd_real_t CutoffRadius;
    fmd_real_t MD_time;
    int time_iteration;
    fmd_real_t delta_t;
    unsigned TotalNoOfParticles;
    unsigned TotalNoOfMolecules;
    fmd_real_t GlobalTemperature;         /* for the "active group" only */
    fmd_bool_t Is_MD_process;
    fmd_bool_t Is_MD_comm_root;
    int LOP_iteration;                    // must be initialized with zero
    int LOP_period;
    MPI_Comm MD_comm;
    fmd_real_t TotalKineticEnergy;
    fmd_real_t TotalPotentialEnergy;
    fmd_real_t TotalMDEnergy;
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
    fmd_bool_t CompLocOrdParam;           // compute local order parameter?
    fmd_SaveConfigMode_t SaveConfigMode;
    FILE *ConfigFilep;
    fmd_real_t WallTimeOrigin;
    int ActiveGroup;
    unsigned ActiveGroupParticlesNum;
    fmd_rtuple_t TotalMomentum;
    fmd_bool_t ParticlesDistributed;
    fmd_bool_t MPI_initialized_by_me;
    h5_dataspaces_t h5_dataspaces;
    int cell_increment;
    int _OldNumberOfParticles;
    int _FileIndex;
    fmd_real_t _OldTotalMDEnergy;
    fmd_real_t _PrevFailedMDEnergy;
};

typedef struct _fmd fmd_t;

/* Functions */

void fmd_box_createGrid(fmd_t *md, fmd_real_t cutoff);
void fmd_dync_setBerendsenThermostatParameter(fmd_t *md, fmd_real_t parameter);
void compLocOrdParam(fmd_t *md);
void createCommunicators(fmd_t *md);
void identifyProcess(fmd_t *md);
void loadStateFile(fmd_t *md, cell_t ***global_grid);
void rescaleVelocities(fmd_t *md);
void _fmd_initialize_grid(cell_t ***grid, cellinfo_t *cinfo, unsigned dim1, unsigned dim2, unsigned dim3);
void fmd_dync_setTimeStep(fmd_t *md, fmd_real_t timeStep);

#endif /* BASE_H */
