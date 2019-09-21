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

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <limits.h>
#include <complex.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "potential.h"

// Global macroes and symbolic constants

#define INDEX(ic, nc)   ((ic)[0] + (nc)[0]*((ic)[1] + (nc)[1]*(ic)[2]))

#define INVERSEINDEX(i, nc, ic)                                        \
    ((ic)[0]= (i)%(nc)[0],                                             \
     (ic)[1]=((i)/(nc)[0])%(nc)[1],                                    \
     (ic)[2]=((i)/(nc)[0])/(nc)[1])

#define SQR(x)          ((x)*(x))
#define CUBE(x)         ((x)*(x)*(x))

#define ITERATE(iv, minv, upv)                                         \
    for ( (iv)[0]=(minv)[0]; (iv)[0]<(upv)[0]; (iv)[0]++ )             \
        for ( (iv)[1]=(minv)[1]; (iv)[1]<(upv)[1]; (iv)[1]++ )         \
            for ( (iv)[2]=(minv)[2]; (iv)[2]<(upv)[2]; (iv)[2]++ )

#define SET_jc_IN_DIRECTION(dd)                                        \
    if (md->ns[dd] == 1)                                               \
    {                                                                  \
        if ((kc[dd] == -1) || (kc[dd] == md->nc[dd]))                  \
            if (md->PBC[dd])                                           \
                jc[dd] = (kc[dd] + md->nc[dd]) % md->nc[dd];           \
            else                                                       \
                continue;                                              \
        else                                                           \
            jc[dd] = kc[dd];                                           \
    }                                                                  \
    else                                                               \
        jc[dd] = kc[dd];

#define ROOTPROCESS(numprocs)       ((numprocs) - 1)
#define K_BOLTZMANN                 8.6173303e-5       // (eV / Kelvin)
#define LIGHT_SPEED                 2.9979245800e+06   // (ang / ps)
#define METER_PER_SECOND            1e-2               // (ang / ps)
#define JOULE_PER_METER2            6.2415091259e-02   // (eV / ang^2)
#define JOULE_PER_METER3_KELVIN2    6.2415091259e-12   // (eV / ang^3 Kelvin^2)
#define JOULE_PER_METER3_KELVIN     6.2415091259e-12   // (eV / ang^3 Kelvin)
#define WATT_PER_METER3_KELVIN      6.2415091259e-24   // (eV / ps ang^3 Kelvin)
#define WATT_PER_METER_KELVIN       6.2415091259e-4    // (eV / ps ang Kelvin)
#define MD_MASS_UNIT                9.6485332907e+03   // x unified atomic mass unit
#define PASCAL                      6.2415091259e-12   // (mass_unit / ang ps^2)
#define MD_CHARGE_UNIT              1.2657711566e-10   // x (electrostatic unit of charge (esu) = statcoulomb)
#define E_CHARGE                    3.7946864629e+00   // electron charge in MD electric charge unit
#define E_MASS                      5.6856300621e-08   // electron mass in MD mass unit
#define HBAR                        6.5821195136e-04   // Planck constant devided by 2*pi in MD units (eV ps)
#define GRAM_PER_CM3                6.2415091259e-05   // x (mass_unit / ang^3)
#define NEWTON_PER_METER            6.2415091259e-02   // x (eV / ang^2)
#define MAX_PATH_LENGTH             256

// error codes
#define ERROR_NC_TOO_SMALL                      1
#define ERROR_UNEXPECTED_PARTICLE_POSITION      2
#define ERROR_UNABLE_OPEN_FILE                  3
#define ERROR_UNSUITABLE_FILE                   4

// typedefs & structs

typedef struct
{
    double x[3];
    double v[3];
    double x_bak[3];
    double v_bak[3];
    float LocOrdParam;
    float x_avgd[3];
    unsigned AtomID;
    unsigned atomkind;
    int GroupID;
    unsigned molkind;
    unsigned MolID;
    unsigned AtomID_local;
} TParticle;

typedef struct TParticleListItem
{
    TParticle P;
    double F[3];
    double FembPrime;
    float LocOrdParamAvg;
    struct TParticleListItem *next_p;
    list_t *neighbors;  // each data pointer in this list points to a mol_atom_neighbor_t
} TParticleListItem;

typedef struct
{
    TParticleListItem *atom;
    unsigned LocalID;
    bondkind_t *bond;
} mol_atom_neighbor_t;

typedef TParticleListItem *TCell;

typedef struct
{
    float x[3];
    float var;
    unsigned atomkind;
} TXYZ_Struct;

typedef struct
{
    double x[3];
    unsigned atomkind;
    int GroupID;
} TPosition_Struct;

typedef struct
{
    TCell ***grid;              // where the particles lie
    int myrank;                 // rank of the local process in MD_comm
    int numprocs;               // number of processes in MD_comm
    int is[3];                  // position of subdomain in the subdomain grid
    int rank_of_lower_subd[3];  // rank of the neighbor processes
    int rank_of_upper_subd[3];
    int ic_start[3];            // width of margin, corresponds to the first
                                // local index in the interior of the subdomain
    int ic_stop[3];             // first local index in the upper margin
    int cell_num[3];            // number of cells in subdomain, including margin
    unsigned cell_num_nonmarg[3];
    int ic_global_firstcell[3]; // global index of the first cell of the subdomain
    unsigned NumberOfParticles;
} TSubDomain;

typedef enum
{
    FMD_SCM_XYZ_PARTICLESNUM,
    FMD_SCM_XYZ_SEPARATE,
    FMD_SCM_CSV,
    FMD_SCM_VTF
} fmd_SaveConfigMode_t;

typedef struct _fmd fmd_t;
typedef struct _fmd_timer fmd_timer_t;

typedef enum
{
    FMD_EVENT_TIMERTICK
} fmd_event_t;

typedef void (*fmd_EventHandler_t)(fmd_t *md, fmd_event_t event, unsigned param);

struct _fmd
{
    TSubDomain SubDomain;
    potsys_t potsys;
    TCell ***global_grid;
    fmd_EventHandler_t eventHandler;
    unsigned timers_num;
    fmd_timer_t *timers;
    fmd_bool_t GlobalGridExists;
    fmd_bool_t boxSizeDetermined;
    fmd_bool_t PBCdetermined;
    double cutoffRadius;
    double mdTime;
    double delta_t;
    unsigned TotalNoOfParticles;
    unsigned TotalNoOfMolecules;
    double GlobalTemperature;
    fmd_bool_t isMDprocess;
    fmd_bool_t Is_MD_comm_root;
    int LOPiteration;               // must be initialized with zero
    MPI_Comm MD_comm;
    double totalKineticEnergy;
    double totalPotentialEnergy;
    double totalMDEnergy;
    int world_rank;
    int world_numprocs;
    double DesiredTemperature;
    int PBC[3];
    int ns[3];                      // number of subdomains = ns[0] x ns[1] x ns[2]
    double l[3];                    // size of the simulation box
    int nc[3];                      // number of grid cells in the simulation box
    double cellh[3];                // size of one single grid cell
    fmd_bool_t UseAutoStep;
    double autoStepSensitivity;
    char saveDirectory[MAX_PATH_LENGTH];
    double BerendsenThermostatParam;
    fmd_bool_t CompLocOrdParam;           // compute local order parameter?
    int locOrdParamPeriod;
    fmd_SaveConfigMode_t saveConfigMode;
    FILE *configFilep;
    double wallTimeOrigin;
    int activeGroup;
    int activeGroupParticlesNum;
    double totalMomentum[3];
    fmd_bool_t ParticlesDistributed;
    fmd_bool_t MPI_initialized_by_me;
    int _oldNumberOfParticles;
    int _fileIndex;
    double _oldTotalMDEnergy;
    double _prevFailedMDEnergy;
};

// Functions

void fmd_subd_init(fmd_t *md);
void fmd_box_createGrid(fmd_t *md, double cutoff);
void fmd_dync_setBerendsenThermostatParameter(fmd_t *md, double parameter);
void cleanGridSegment(TCell ***grid, int ic_from[3], int ic_to[3]);
void compLocOrdParam(fmd_t *md);
void createCommunicators(fmd_t *md);
TCell ***createGrid(int cell_num[3]);
int getListLength(TParticleListItem *root_p);
void identifyProcess(fmd_t *md);
void handleFileOpenError(FILE *fp, char *filename);
void loadStateFile(fmd_t *md, TCell ***global_grid);
void rescaleVelocities(fmd_t *md);
void restoreBackups(fmd_t *md);
void fmd_insertInList(TParticleListItem **root_pp, TParticleListItem *item_p);
void removeFromList(TParticleListItem **item_pp);

//

extern const int fmd_ThreeZeros[3];

#endif /* BASE_H */
