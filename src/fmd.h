/*
  fmd.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

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

#ifndef FMD_H
#define FMD_H

#define FMD_TRUE 1
#define FMD_FALSE 0

// typedefs

typedef struct _fmd fmd_t;

typedef struct _fmd_pot fmd_pot_t;

typedef enum
{
    FMD_SCM_XYZ_PARTICLESNUM,
    FMD_SCM_XYZ_SEPARATE,
    FMD_SCM_CSV,
    FMD_SCM_VTF
} fmd_SaveConfigMode_t;

typedef enum
{
    FMD_BOND_HARMONIC
} fmd_bond_t;

typedef enum
{
    FMD_TURI_CUSTOM,
    FMD_TURI_TTM
} fmd_turi_t;

typedef char *fmd_string_t;
typedef int fmd_bool_t;

typedef enum
{
    FMD_EVENT_TIMERTICK
} fmd_event_t;

typedef void (*fmd_EventHandler_t)(fmd_t *md, fmd_event_t event, unsigned param);

// functions

unsigned fmd_bond_addKind(fmd_t *md, fmd_bond_t cat, double coeffs[]);
void fmd_bond_apply(fmd_t *md, unsigned bondkind, unsigned molkind,
  unsigned atom1, unsigned atom2);

unsigned fmd_molecule_addKind(fmd_t *md, fmd_string_t name, unsigned AtomsNum,
  unsigned AtomKinds[], double AtomPositions[][3]);

void fmd_matt_addVelocity(fmd_t *md, int GroupID, double vx, double vy, double vz);
void fmd_matt_setActiveGroup(fmd_t *md, int GroupID);
void fmd_matt_setDesiredTemperature(fmd_t *md, double DesiredTemperature);
void fmd_matt_makeCuboidSC(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, unsigned atomkind, int GroupID);
void fmd_matt_makeCuboidSC_alloy(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, double *proportions, int GroupID);
void fmd_matt_makeCuboidBCC(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, unsigned atomkind, int GroupID);
void fmd_matt_makeCuboidBCC_alloy(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, double *proportions, int GroupID);
void fmd_matt_makeCuboidFCC(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, unsigned atomkind, int GroupID);
void fmd_matt_makeCuboidFCC_alloy(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double LatticeParameter, double *proportions, int GroupID);
void fmd_matt_scatterMolecule(fmd_t *md, unsigned molkind, double xa,
  double ya, double za, double xb, double yb, double zb, unsigned num,
  int GroupID);
void fmd_matt_saveConfiguration(fmd_t *md);
double fmd_matt_getTotalEnergy(fmd_t *md);
double fmd_matt_getGlobalTemperature(fmd_t *md);
void fmd_matt_distribute(fmd_t *md);
void fmd_matt_giveTemperature(fmd_t *md, int GroupID);
void fmd_matt_setAtomKinds(fmd_t *md, unsigned number, const fmd_string_t names[], const double masses[]);

void fmd_box_setPBC(fmd_t *md, fmd_bool_t PBCx, fmd_bool_t PBCy, fmd_bool_t PBCz);
void fmd_box_setSize(fmd_t *md, double sx, double sy, double sz);
void fmd_box_setSubDomains(fmd_t *md, int dimx, int dimy, int dimz);
void fmd_box_createGrid(fmd_t *md, double cutoff);

void fmd_io_setSaveDirectory(fmd_t *md, fmd_string_t directory);
void fmd_io_setSaveConfigMode(fmd_t *md, fmd_SaveConfigMode_t mode);
void fmd_io_printf(fmd_t *md, const fmd_string_t restrict format, ...);
void fmd_io_loadState(fmd_t *md, fmd_string_t filepath, fmd_bool_t useTime);
void fmd_io_saveState(fmd_t *md, fmd_string_t filename);

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_t *md, fmd_string_t filePath);
double fmd_pot_eam_getLatticeParameter(fmd_t *md, fmd_pot_t *pot, fmd_string_t element);
double fmd_pot_eam_getCutoffRadius(fmd_t *md, fmd_pot_t *pot);
void fmd_pot_setCutoffRadius(fmd_t *md, double cutoff);
fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  double sigma, double epsilon, double cutoff);
fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  double D0, double alpha, double r0, double cutoff);
void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);

void fmd_subd_init(fmd_t *md);
void fmd_subd_free(fmd_t *md);

unsigned fmd_timer_makeSimple(fmd_t *md, double start, double interval, double stop);

double fmd_proc_getWallTime(fmd_t *md);
fmd_bool_t fmd_proc_isMD(fmd_t *md);
fmd_bool_t fmd_proc_isRoot(fmd_t *md);

fmd_t *fmd_create();
void fmd_free(fmd_t *md);
void fmd_setEventHandler(fmd_t *md, fmd_EventHandler_t func);

double fmd_dync_getTimeStep(fmd_t *md);
void fmd_dync_setTimeStep(fmd_t *md, double TimeStep);
double fmd_dync_getTime(fmd_t *md);
void fmd_dync_updateForces(fmd_t *md);
void fmd_dync_updateForcesLJ(fmd_t *md);
void fmd_dync_incTime(fmd_t *md);
void fmd_dync_setBerendsenThermostatParameter(fmd_t *md, double parameter);
void fmd_dync_VelocityVerlet_startStep(fmd_t *md, fmd_bool_t UseThermostat);
int fmd_dync_VelocityVerlet_finishStep(fmd_t *md);
void fmd_dync_equilibrate(fmd_t *md, int GroupID, double duration,
  double timestep, double strength, double temperature);

unsigned fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz);

#endif /* FMD_H */
