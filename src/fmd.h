/*
  fmd.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, and other contributors

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

/* typedefs */

typedef enum
{
    FMD_SCM_XYZ_PARTICLESNUM,
    FMD_SCM_XYZ_SEPARATE,
    FMD_SCM_CSV,
    FMD_SCM_VTF
} fmd_SaveConfigMode_t;

typedef enum
{
    FMD_EVENT_TIMER_TICK,
    FMD_EVENT_FIELD_UPDATE
} fmd_event_t;

typedef enum
{
    FMD_BOND_HARMONIC
} fmd_bond_t;

typedef enum
{
    FMD_TURI_CUSTOM,
    FMD_TURI_TTM_TYPE1
} fmd_turi_t;

typedef enum
{
    FMD_FIELD_MASS,
    FMD_FIELD_VCM,
    FMD_FIELD_TEMPERATURE,
    FMD_FIELD_NUMBER,
    FMD_FIELD_NUMBER_DENSITY
} fmd_field_t;

typedef enum
{
    FMD_TTM_TIMESTEP_RATIO_CONSTANT
} fmd_ttm_timestep_ratio_t;

typedef enum
{
    FMD_TTM_TE_CONSTANT
} fmd_ttm_Te_t;

typedef struct _fmd fmd_t;

typedef struct _fmd_pot fmd_pot_t;

typedef double fmd_real_t;

typedef fmd_real_t fmd_rtuple_t[3];

typedef int fmd_handle_t;

typedef char *fmd_string_t;

typedef int fmd_bool_t;

typedef struct _fmd_params fmd_params_t;

typedef void (*fmd_EventHandler_t)(fmd_t *md, fmd_event_t event, fmd_params_t *params);

typedef struct
{
    fmd_handle_t timer;
} fmd_event_params_timer_tick_t;

typedef struct
{
    fmd_handle_t turi;
    fmd_handle_t field;
} fmd_event_params_field_update_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_params_Te_constant_t;

typedef struct
{
    fmd_real_t gamma;
} fmd_ttm_params_heat_capacity_linear_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_params_heat_conductivity_constant_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_params_coupling_factor_constant_t;

typedef struct
{
    unsigned value;
} fmd_ttm_params_timestep_ratio_constant_t;

/* functions */

fmd_handle_t fmd_bond_addKind(fmd_t *md, fmd_bond_t cat, fmd_real_t coeffs[]);
void fmd_bond_apply(fmd_t *md, fmd_handle_t bondkind, fmd_handle_t molkind,
  unsigned atom1, unsigned atom2);

fmd_handle_t fmd_molecule_addKind(fmd_t *md, fmd_string_t name, unsigned AtomsNum,
  unsigned AtomKinds[], fmd_rtuple_t AtomPositions[]);

void fmd_matt_addVelocity(fmd_t *md, int GroupID, fmd_real_t vx, fmd_real_t vy, fmd_real_t vz);
void fmd_matt_setDesiredTemperature(fmd_t *md, fmd_real_t DesiredTemperature);
void fmd_matt_makeCuboidSC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID);
void fmd_matt_makeCuboidSC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID);
void fmd_matt_makeCuboidBCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID);
void fmd_matt_makeCuboidBCC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID);
void fmd_matt_makeCuboidFCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, unsigned atomkind, int GroupID);
void fmd_matt_makeCuboidFCC_alloy(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t LatticeParameter, fmd_real_t *proportions, int GroupID);
void fmd_matt_scatterMolecule(fmd_t *md, fmd_handle_t molkind, fmd_real_t xa,
  fmd_real_t ya, fmd_real_t za, fmd_real_t xb, fmd_real_t yb, fmd_real_t zb, unsigned num,
  int GroupID);
void fmd_matt_saveConfiguration(fmd_t *md);
fmd_real_t fmd_matt_getTotalEnergy(fmd_t *md);
fmd_real_t fmd_matt_getGroupTemperature(fmd_t *md);
void fmd_matt_giveTemperature(fmd_t *md, int GroupID);
void fmd_matt_setAtomKinds(fmd_t *md, unsigned number, const fmd_string_t names[], const fmd_real_t masses[]);

void fmd_box_setPBC(fmd_t *md, fmd_bool_t PBCx, fmd_bool_t PBCy, fmd_bool_t PBCz);
void fmd_box_setSize(fmd_t *md, fmd_real_t sx, fmd_real_t sy, fmd_real_t sz);
void fmd_box_setSubDomains(fmd_t *md, int dimx, int dimy, int dimz);
void fmd_box_createGrid(fmd_t *md, fmd_real_t cutoff);

void fmd_io_setSaveDirectory(fmd_t *md, fmd_string_t directory);
void fmd_io_setSaveConfigMode(fmd_t *md, fmd_SaveConfigMode_t mode);
void fmd_io_printf(fmd_t *md, const fmd_string_t restrict format, ...);
void fmd_io_loadState(fmd_t *md, fmd_string_t filepath, fmd_bool_t useTime);
void fmd_io_saveState(fmd_t *md, fmd_string_t filename);

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_t *md, fmd_string_t filePath);
fmd_real_t fmd_pot_eam_getLatticeParameter(fmd_t *md, fmd_pot_t *pot, fmd_string_t element);
fmd_real_t fmd_pot_eam_getCutoffRadius(fmd_t *md, fmd_pot_t *pot);
void fmd_pot_setCutoffRadius(fmd_t *md, fmd_real_t cutoff);
fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  fmd_real_t sigma, fmd_real_t epsilon, fmd_real_t cutoff);
fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  fmd_real_t D0, fmd_real_t alpha, fmd_real_t r0, fmd_real_t cutoff);
void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);

fmd_handle_t fmd_timer_makeSimple(fmd_t *md, fmd_real_t start, fmd_real_t interval, fmd_real_t stop);

fmd_real_t fmd_proc_getWallTime(fmd_t *md);
fmd_bool_t fmd_proc_isMD(fmd_t *md);
fmd_bool_t fmd_proc_isRoot(fmd_t *md);

fmd_t *fmd_create();
void fmd_free(fmd_t *md);
void fmd_setEventHandler(fmd_t *md, fmd_EventHandler_t func);

fmd_real_t fmd_dync_getTimeStep(fmd_t *md);
fmd_real_t fmd_dync_getTime(fmd_t *md);
void fmd_dync_equilibrate(fmd_t *md, int GroupID, fmd_real_t duration,
  fmd_real_t timestep, fmd_real_t strength, fmd_real_t temperature);
void fmd_dync_integrate(fmd_t *md, int GroupID, fmd_real_t duration, fmd_real_t timestep);

fmd_handle_t fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz, fmd_real_t starttime, fmd_real_t stoptime);
fmd_handle_t fmd_field_add(fmd_t *md, fmd_handle_t turi, fmd_field_t cat, fmd_real_t interval);
void fmd_field_save_as_hdf5(fmd_t *md, fmd_handle_t turi, fmd_handle_t field, fmd_string_t path);

void fmd_ttm_setHeatCapacity(fmd_t *md, fmd_handle_t turi, fmd_params_t *params);
void fmd_ttm_setHeatConductivity(fmd_t *md, fmd_handle_t turi, fmd_params_t *params);
void fmd_ttm_setCouplingFactor(fmd_t *md, fmd_handle_t turi, fmd_params_t *params);
void fmd_ttm_setElectronTemperature(fmd_t *md, fmd_handle_t turi, fmd_ttm_Te_t cat, fmd_params_t *params);
void fmd_ttm_setTimestepRatio(fmd_t *md, fmd_handle_t turi, fmd_ttm_timestep_ratio_t cat, fmd_params_t *params);

#endif /* FMD_H */
