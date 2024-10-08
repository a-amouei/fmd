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

#include <stdbool.h>

#define FMD_GROUP_ALL -1

/* typedefs */

typedef void *fmd_pointer_t;

typedef enum
{
    FMD_SCM_XYZ_ATOMSNUM,
    FMD_SCM_XYZ_SEPARATE,
    FMD_SCM_CSV,
    FMD_SCM_VTF
} fmd_SaveConfigMode_t;

typedef enum
{
    FMD_EVENT_TIMER_TICK,
    FMD_EVENT_FIELD_UPDATE,
    FMD_EVENT_ERROR
} fmd_event_t;

typedef enum
{
    FMD_TURI_CUSTOM,
    FMD_TURI_TTM_TYPE1,
    FMD_TURI_TTM_TYPE2
} fmd_turi_t;

typedef enum
{
    FMD_FIELD_MASS,
    FMD_FIELD_VCM,
    FMD_FIELD_TEMPERATURE,
    FMD_FIELD_NUMBER,
    FMD_FIELD_NUMBER_DENSITY,
    FMD_FIELD_TTM_TE,
    FMD_FIELD_TTM_XI
} fmd_field_t;

typedef struct _fmd fmd_t;

typedef struct _fmd_array3s fmd_array3s_t;

typedef struct _fmd_pot fmd_pot_t;

typedef double fmd_real_t;

typedef fmd_real_t fmd_rtriple_t[3];

typedef unsigned fmd_utriple_t[3];

typedef void ***fmd_array3_t;

typedef int fmd_handle_t;

typedef char *fmd_string_t;

typedef struct _fmd_params fmd_params_t;

typedef void (*fmd_EventHandler_t)(fmd_t *md, fmd_event_t event, void *usp, fmd_params_t *params);

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

typedef struct
{
    fmd_error_t error;
    bool major;
    fmd_string_t source;
    fmd_string_t func;
    int line;
    fmd_pointer_t p1;
    fmd_pointer_t p2;
    fmd_pointer_t p3;
} fmd_event_params_error_t;

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
    fmd_real_t gamma;
} fmd_ttm_heat_capacity_linear_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_heat_conductivity_constant_t;

typedef struct
{
    fmd_real_t v;
    fmd_real_t A;
    fmd_real_t B;
} fmd_ttm_heat_conductivity_zhigilei_t;

typedef struct
{
    fmd_real_t value;
} fmd_ttm_coupling_factor_constant_t;

typedef struct
{
    fmd_real_t fluence;
    fmd_real_t reflectance;
    fmd_real_t t0;
    fmd_real_t duration;
    fmd_real_t AbsorptionDepth;
} fmd_ttm_laser_gaussian_t;

/* functions */

void fmd_matt_addVelocity(fmd_t *md, int GroupID, fmd_real_t vx, fmd_real_t vy, fmd_real_t vz);
void fmd_matt_translate(fmd_t *md, int GroupID, fmd_real_t dx, fmd_real_t dy, fmd_real_t dz);
void fmd_matt_makeCuboidSC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, unsigned atomkind, int GroupID, fmd_real_t temp);
void fmd_matt_makeCuboidSC_mix(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp);
void fmd_matt_makeCuboidBCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, unsigned atomkind, int GroupID, fmd_real_t temp);
void fmd_matt_makeCuboidBCC_mix(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp);
void fmd_matt_makeCuboidFCC(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, unsigned atomkind, int GroupID, fmd_real_t temp);
void fmd_matt_makeCuboidFCC_mix(fmd_t *md, fmd_real_t x, fmd_real_t y, fmd_real_t z,
  int dimx, int dimy, int dimz, fmd_real_t lp, fmd_real_t ratio[], int GroupID, fmd_real_t temp);
void fmd_matt_saveConfiguration(fmd_t *md);
fmd_real_t fmd_matt_getKineticEnergy(fmd_t *md);
fmd_real_t fmd_matt_getTotalEnergy(fmd_t *md);
fmd_real_t fmd_matt_getTemperature(fmd_t *md);
void fmd_matt_getMomentum(fmd_t *md, fmd_rtriple_t out);
void fmd_matt_giveMaxwellDistribution(fmd_t *md, int GroupID, fmd_real_t temp);
void fmd_matt_setAtomKinds(fmd_t *md, unsigned number, const fmd_string_t names[], const fmd_real_t masses[]);
void fmd_matt_findLimits(fmd_t *md, int GroupID, fmd_rtriple_t LowerLimit, fmd_rtriple_t UpperLimit);
void fmd_matt_changeGroupID(fmd_t *md, int old, int new);

void fmd_box_setPBC(fmd_t *md, bool PBCx, bool PBCy, bool PBCz);
void fmd_box_setSize(fmd_t *md, fmd_real_t sx, fmd_real_t sy, fmd_real_t sz);
bool fmd_box_setSubdomains(fmd_t *md, int dimx, int dimy, int dimz);

void fmd_io_setSaveDirectory(fmd_t *md, fmd_string_t directory);
void fmd_io_setSaveConfigMode(fmd_t *md, fmd_SaveConfigMode_t mode);
void fmd_io_printf(fmd_t *md, const fmd_string_t restrict format, ...);
void fmd_io_loadState(fmd_t *md, fmd_string_t path, bool UseTime);
void fmd_io_saveState(fmd_t *md, fmd_string_t filename);
void fmd_io_setShowErrorMessages(fmd_t *md, bool show);

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_t *md, fmd_string_t path);
fmd_real_t fmd_pot_eam_getLatticeParameter(fmd_t *md, fmd_pot_t *pot, fmd_string_t element);
fmd_real_t fmd_pot_eam_getCutoffRadius(fmd_t *md, fmd_pot_t *pot);
fmd_pot_t *fmd_pot_lj_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  fmd_real_t sigma, fmd_real_t epsilon, fmd_real_t cutoff);
fmd_pot_t *fmd_pot_morse_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2,
  fmd_real_t D0, fmd_real_t alpha, fmd_real_t r0, fmd_real_t cutoff);
void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot);

fmd_handle_t fmd_timer_makeSimple(fmd_t *md, fmd_real_t start, fmd_real_t interval, fmd_real_t stop);

fmd_real_t fmd_proc_getWallTime(fmd_t *md);
bool fmd_proc_hasSubdomain(fmd_t *md);
bool fmd_proc_isRoot(fmd_t *md);
void fmd_proc_setNumThreads(fmd_t *md, int num);
int fmd_proc_getCellIncrement(fmd_t *md);
void fmd_proc_setCellIncrement(fmd_t *md , int incr);

fmd_t *fmd_create();
void fmd_free(fmd_t *md);
void fmd_setEventHandler(fmd_t *md, void *usp, fmd_EventHandler_t func);

fmd_real_t fmd_dync_getTimestep(fmd_t *md);
fmd_real_t fmd_dync_getTime(fmd_t *md);
void fmd_dync_equilibrate(fmd_t *md, int GroupID, fmd_real_t duration,
  fmd_real_t timestep, fmd_real_t tau, fmd_real_t temperature);
void fmd_dync_integrate(fmd_t *md, int GroupID, fmd_real_t duration, fmd_real_t timestep);

fmd_handle_t fmd_turi_add(fmd_t *md, fmd_turi_t cat, int dimx, int dimy, int dimz, fmd_real_t starttime, fmd_real_t stoptime);
fmd_handle_t fmd_field_find(fmd_t *md, fmd_handle_t turi, fmd_field_t cat);
fmd_handle_t fmd_field_add(fmd_t *md, fmd_handle_t turi, fmd_field_t cat, fmd_real_t interval);
void fmd_field_save_as_hdf5(fmd_t *md, fmd_handle_t turi, fmd_handle_t field, fmd_string_t filename);
fmd_array3s_t *fmd_field_getArray(fmd_t *md, fmd_handle_t turi, fmd_handle_t field,
  fmd_array3_t *array, fmd_utriple_t dims);

void _fmd_ttm_setHeatCapacity_linear(fmd_t *md, fmd_handle_t turi, fmd_ttm_heat_capacity_linear_t c);
void _fmd_ttm_setHeatCapacity_file(fmd_t *md, fmd_handle_t turi, fmd_string_t path);
void _fmd_ttm_setHeatConductivity_constant1(fmd_t *md, fmd_handle_t turi, fmd_real_t k);
void _fmd_ttm_setHeatConductivity_constant2(fmd_t *md, fmd_handle_t turi, fmd_ttm_heat_conductivity_constant_t k);
void _fmd_ttm_setHeatConductivity_zhigilei(fmd_t *md, fmd_handle_t turi, fmd_ttm_heat_conductivity_zhigilei_t k);
void _fmd_ttm_setCouplingFactor_constant1(fmd_t *md, fmd_handle_t turi, fmd_real_t g);
void _fmd_ttm_setCouplingFactor_constant2(fmd_t *md, fmd_handle_t turi, fmd_ttm_coupling_factor_constant_t g);
void _fmd_ttm_setCouplingFactor_file(fmd_t *md, fmd_handle_t turi, fmd_string_t path);
void fmd_ttm_setElectronTemperature(fmd_t *md, fmd_handle_t turi, fmd_real_t Te);
void fmd_ttm_setTimestepRatio(fmd_t *md, fmd_handle_t turi, int ratio);
void fmd_ttm_setCellActivationFraction(fmd_t *md, fmd_handle_t turi, fmd_real_t value);
void _fmd_ttm_setLaserSource_gaussian(fmd_t *md, fmd_handle_t turi, fmd_ttm_laser_gaussian_t laser);

#define fmd_ttm_setHeatCapacity(md, turi, c) \
  _Generic((c), fmd_ttm_heat_capacity_linear_t: _fmd_ttm_setHeatCapacity_linear, \
                fmd_string_t: _fmd_ttm_setHeatCapacity_file)(md, turi, c)

#define fmd_ttm_setHeatConductivity(md, turi, k) \
  _Generic((k), fmd_ttm_heat_conductivity_constant_t: _fmd_ttm_setHeatConductivity_constant2, \
                fmd_real_t: _fmd_ttm_setHeatConductivity_constant1, \
                fmd_ttm_heat_conductivity_zhigilei_t: _fmd_ttm_setHeatConductivity_zhigilei)(md, turi, k)

#define fmd_ttm_setCouplingFactor(md, turi, g) \
  _Generic((g), fmd_ttm_coupling_factor_constant_t: _fmd_ttm_setCouplingFactor_constant2, \
                fmd_real_t: _fmd_ttm_setCouplingFactor_constant1, \
                fmd_string_t: _fmd_ttm_setCouplingFactor_file)(md, turi, g)

#define fmd_ttm_setLaserSource(md, turi, laser) \
  _Generic((laser), fmd_ttm_laser_gaussian_t: _fmd_ttm_setLaserSource_gaussian)(md, turi, laser)

void fmd_array3s_free(fmd_array3s_t *array);

int fmd_version_getMajor();
int fmd_version_getMinor();
int fmd_version_getRevision();
char fmd_version_getType();
fmd_string_t fmd_version_getString();

#endif /* FMD_H */
