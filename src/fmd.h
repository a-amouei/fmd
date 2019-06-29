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

// types

typedef struct fmd_t fmd_t;

typedef struct fmd_pot_t fmd_pot_t;

typedef enum
{
    SCM_XYZ_PARTICLESNUM,
    SCM_XYZ_SEPARATE,
    SCM_CSV,
    SCM_VTF
} fmd_SaveConfigMode_t;

typedef char *fmd_string_t;

// functions

void fmd_matt_addVelocity(fmd_t *md, int groupID, double vx, double vy, double vz);
void fmd_matt_setActiveGroup(fmd_t *md, int groupID);
void fmd_matt_setDesiredTemperature(fmd_t *md, double DesiredTemperature);
void fmd_matt_makeCuboidFCC(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double latticeParameter, int elementID, int groupID);
void fmd_matt_makeCuboidFCC_alloy(fmd_t *md, double x, double y, double z,
  int dimx, int dimy, int dimz, double latticeParameter, double *proportions, int groupID);
void fmd_matt_saveConfiguration(fmd_t *md);
double fmd_matt_getTotalEnergy(fmd_t *md);
double fmd_matt_getGlobalTemperature(fmd_t *md);
void fmd_matt_distribute(fmd_t *md);
void fmd_matt_giveTemperature(fmd_t *md, int groupID);

void fmd_box_setPBC(fmd_t *md, int PBCx, int PBCy, int PBCz);
void fmd_box_setSize(fmd_t *md, double sx, double sy, double sz);
void fmd_box_setSubDomains(fmd_t *md, int dimx, int dimy, int dimz);
void fmd_box_createGrid(fmd_t *md, double cutoff);

void fmd_io_setSaveDirectory(fmd_t *md, char *directory);
void fmd_io_setSaveConfigMode(fmd_t *md, fmd_SaveConfigMode_t mode);
void fmd_io_printf(fmd_t *md, const char * restrict format, ...);
void fmd_io_loadState(fmd_t *md, char *filepath, int useTime);
void fmd_io_saveState(fmd_t *md, char *filename);

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_t *md, char *filePath);
double fmd_pot_eam_getLatticeParameter(fmd_t *md, fmd_pot_t *pot, fmd_string_t element);
double fmd_pot_eam_getCutoffRadius(fmd_t *md, fmd_pot_t *pot);
void fmd_pot_setCutoffRadius(fmd_t *md, double cutoff);
void fmd_pot_setAtomKinds(fmd_t *md, unsigned number, const fmd_string_t names[], const double masses[]);
fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  double sigma, double epsilon, double cutoff);
fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
                               double D0, double alpha, double r0, double cutoff);
void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);

void fmd_subd_init(fmd_t *md);
void fmd_subd_free(fmd_t *md);

double fmd_proc_getWallTime(fmd_t *md);
int fmd_proc_isMD(fmd_t *md);

fmd_t *fmd_sys_create();
void fmd_sys_free(fmd_t *md, int finalizeMPI);

double fmd_dync_getTimeStep(fmd_t *md);
void fmd_dync_setTimeStep(fmd_t *md, double timeStep);
double fmd_dync_getTime(fmd_t *md);
void fmd_dync_updateForces(fmd_t *md);
void fmd_dync_updateForcesLJ(fmd_t *md);
void fmd_dync_incTime(fmd_t *md);
void fmd_dync_setBerendsenThermostatParameter(fmd_t *md, double parameter);
void fmd_dync_velocityVerlet_takeFirstStep(fmd_t *md, int useThermostat);
int fmd_dync_velocityVerlet_takeLastStep(fmd_t *md);
void fmd_dync_equilibrate(fmd_t *md, int groupID, double duration, double strength);

#endif /* FMD_H */
