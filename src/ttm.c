/*
  ttm.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2022 Arham Amouye Foumani

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
#include <tgmath.h>
#include "ttm.h"
#include "fmd-private.h"
#include "matter.h"
#include "general.h"
#include "turi.h"
#include "turi-ghost.h"
#include "types.h"
#include "cspline.h"

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

void _fmd_ttm_getReady(fmd_t *md)
{
    turi_t *t = md->active_ttm_turi;
    ttm_t *ttm = t->ttm;

    fmd_rtuple_t ll, ul;

    fmd_matt_findLimits(md, md->ActiveGroup, ll, ul);
    ttm->frontsurf = ll[DIM-1];

    fmd_real_t vol = 1.0;

    for (int d=0; d<DIM; d++)
        vol *= ul[d] - ll[d];

    ttm->initial_atoms_num = md->GroupParticlesNum / vol * t->tcell_volume;
    ttm->min_atoms_num = ttm->initial_atoms_num * ttm->CellActivFrac;
}

static void type2_1d_pack(fmd_t *md, turi_t *t, bool SendToUp, fmd_pointer_t sendbuf, int *pos)
{
    ttm_t *ttm = t->ttm;
    *pos = 0;

    int i = SendToUp ? t->itc_stop[2] - 1 : t->itc_start_owned[2];

    MPI_Pack(&ttm->num_1d[i], 1, MPI_UNSIGNED, sendbuf, INT_MAX, pos, md->MD_comm);
    MPI_Pack(&ttm->Te_1d[i],  1, FMD_MPI_REAL, sendbuf, INT_MAX, pos, md->MD_comm);
    MPI_Pack(&ttm->Ke_1d[i],  1, FMD_MPI_REAL, sendbuf, INT_MAX, pos, md->MD_comm);
}

static void type1_1d_pack(fmd_t *md, turi_t *t, bool SendToUp, fmd_pointer_t sendbuf, int *pos)
{
    ttm_t *ttm = t->ttm;
    *pos = 0;

    int i = SendToUp ? t->itc_stop[2] - 1 : t->itc_start_owned[2];

    MPI_Pack(&ttm->num_1d[i], 1, MPI_UNSIGNED, sendbuf, INT_MAX, pos, md->MD_comm);
    MPI_Pack(&ttm->Te_1d[i],  1, FMD_MPI_REAL, sendbuf, INT_MAX, pos, md->MD_comm);
}

static void type2_1d_unpack(fmd_t *md, turi_t *t, bool SendToUp, fmd_pointer_t recvbuf)
{
    ttm_t *ttm = t->ttm;
    int pos = 0;

    int i = SendToUp ? 0 : t->itc_stop[2];

    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->num_1d[i], 1, MPI_UNSIGNED, md->MD_comm);
    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->Te_1d[i],  1, FMD_MPI_REAL, md->MD_comm);
    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->Ke_1d[i],  1, FMD_MPI_REAL, md->MD_comm);
}

static void type1_1d_unpack(fmd_t *md, turi_t *t, bool SendToUp, fmd_pointer_t recvbuf)
{
    ttm_t *ttm = t->ttm;
    int pos = 0;

    int i = SendToUp ? 0 : t->itc_stop[2];

    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->num_1d[i], 1, MPI_UNSIGNED, md->MD_comm);
    MPI_Unpack(recvbuf, INT_MAX, &pos, &ttm->Te_1d[i],  1, FMD_MPI_REAL, md->MD_comm);
}

static size_t type1_1d_calc_tghost_buffsize(fmd_t *md)
{
    size_t size;
    int i;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &i);
    size = i;
    MPI_Pack_size(1, FMD_MPI_REAL, md->MD_comm, &i);
    size += i;

    return size;
}

static size_t type2_1d_calc_tghost_buffsize(fmd_t *md)
{
    size_t size;
    int i;

    MPI_Pack_size(1, MPI_UNSIGNED, md->MD_comm, &i);
    size = i;
    MPI_Pack_size(1, FMD_MPI_REAL, md->MD_comm, &i);
    size += 2*i;

    return size;
}

static void type2_update_Ce_Ke_G(turi_t *t, ttm_t *ttm)
{
    for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
    {
        int khi, klo, k;
        double a, b, h;

        if (ttm->Te_1d[i] < 0.) continue; /* nothing to do with deactivated cells */

        assert(ttm->Te_1d[i] > ttm->T_C[0]);
        assert(ttm->Te_1d[i] < ttm->T_C[ttm->nC-1]);

        FIND_klo_AND_ETC(ttm->T_C, ttm->Te_1d[i], ttm->nC);

#ifdef USE_CSPLINE
        ttm->Ce_1d[i] = SPLINE_VAL(a, b, ttm->Ct, klo, khi, ttm->C_DD, h);
#else
        ttm->Ce_1d[i] = a * ttm->Ct[klo] + b * ttm->Ct[khi];
#endif

        ttm->Ke_1d[i] = 1./3. * ttm->Ce_1d[i] * ttm->Kzh.v2 /
                        (ttm->Kzh.A * sqrr(ttm->Te_1d[i]) + ttm->Kzh.B * ttm->Ti_1d[i]);

        assert(ttm->Te_1d[i] > ttm->T_G[0]);
        assert(ttm->Te_1d[i] < ttm->T_G[ttm->nG-1]);

        FIND_klo_AND_ETC(ttm->T_G, ttm->Te_1d[i], ttm->nG);

#ifdef USE_CSPLINE
        ttm->G_1d[i] = SPLINE_VAL(a, b, ttm->Gt, klo, khi, ttm->G_DD, h);
#else
        ttm->G_1d[i] = a * ttm->Gt[klo] + b * ttm->Gt[khi];
#endif
    }
}

static void presolve_common_types_1_2(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    if (t->has_upper_lower_owner_procs[2])
        _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp);

    /* xi intialization, cell-activation, cell-deactivation */
    for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
    {
        ttm->xi_1d[i] = 0.0;

        /* cell activation */
        if (ttm->Te_1d[i] < 0.0 && ttm->num_1d[i] >= ttm->min_atoms_num)
        {
            int count = 0;
            fmd_real_t sum = 0.0;

            int im1 = (i == t->itc_start_owned[2]) ? 0 : i-1;

            if (ttm->Te_1d[im1] >= 0.0)
            {
                sum += ttm->Te_1d[im1] * ttm->num_1d[im1];
                count += ttm->num_1d[im1];
            }

            if (ttm->Te_1d[i+1] >= 0.0)
            {
                sum += ttm->Te_1d[i+1] * ttm->num_1d[i+1];
                count += ttm->num_1d[i+1];
            }

            ttm->Te_1d[i] = (count > 0) ? sum/count : ttm->Ti_1d[i];
        }

        if (ttm->num_1d[i] < ttm->min_atoms_num)
            ttm->Te_1d[i] = ttm->Te2_1d[i] = -1.0; /* cell deactivation */
    }
}

static void ttm_type2_presolve_1d(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    presolve_common_types_1_2(md, t, ttm);

    type2_update_Ce_Ke_G(t, ttm);

    if (t->has_upper_lower_owner_procs[2]) /* again! since ghost cells may have been activated/deactivated */
        _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp); /* also transfers correct Ke data */
}

static void ttm_type1_presolve_1d(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    presolve_common_types_1_2(md, t, ttm);

    if (t->has_upper_lower_owner_procs[2]) /* again! since ghost cells may have been activated/deactivated */
        _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp);
}

static void ttm_type1_solve_1d(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    bool CalcSource = fabs(md->time - ttm->laser_t0) < ttm->laser_tdiff ? true : false;
    fmd_real_t dt = md->timestep / ttm->timestep_ratio;

    /* time loop */
    for (int j=0; j < ttm->timestep_ratio; j++)
    {
        fmd_real_t source_spcind;

        if (CalcSource)
            source_spcind = ttm->laser_factor_constant *
                            exp(ttm->laser_m_2sig2_inv * sqrr(md->time + j * dt - ttm->laser_t0));

        /* spatial loop */
        for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
        {
            if (ttm->Te_1d[i] < 0) continue;  /* nothing to be done with a deactivated cell */

            fmd_real_t source;

            if (CalcSource)
            {
                fmd_real_t z = (i - t->itc_glob_to_loc[2] + 0.5) * t->tcellh[2];
                source = source_spcind * exp(ttm->laser_m_absdepth_inv * (z - ttm->frontsurf));
            }
            else
                source = 0.0;

            int im1 = (i == t->itc_start_owned[2]) ? 0 : i-1;
            int ilo = (ttm->Te_1d[im1] < 0) ? i : im1;
            int ihi = (ttm->Te_1d[i+1] < 0) ? i : i+1;

            ttm->Te2_1d[i] = ttm->Te_1d[i] + dt/(ttm->C_gamma * ttm->Te_1d[i]) * (
              ttm->K * (ttm->Te_1d[ihi]-2*ttm->Te_1d[i]+ttm->Te_1d[ilo])/ttm->dz2
              - ttm->G * (ttm->Te_1d[i] - ttm->Ti_1d[i]) + source );

            ttm->xi_1d[i] += ttm->Te2_1d[i];
        } /* end of spatial loop */

        fmd_real_t *tempo = ttm->Te2_1d;
        ttm->Te2_1d = ttm->Te_1d;
        ttm->Te_1d = tempo;

        if (t->has_upper_lower_owner_procs[2] && j < ttm->timestep_ratio-1)
            _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp);
    } /* end of time loop */

    if (ttm->Te2_1d != ((fmd_real_t ***)ttm->Te_aux.data)[0][0])
    {
        fmd_array3s_t tempo = ttm->Te_aux;
        ttm->Te_aux = t->fields[ttm->iTe].data;
        t->fields[ttm->iTe].data = tempo;
    }

    for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
    {
        if (ttm->num_1d[i] >= ttm->min_atoms_num)
            ttm->xi_1d[i] = ttm->G * t->tcell_volume *
                            (ttm->xi_1d[i] / ttm->timestep_ratio - ttm->Ti_1d[i]) /
                            (3 * K_BOLTZMANN * ttm->num_1d[i] * ttm->Ti_1d[i]);
    }
}

static void ttm_type2_solve_1d(fmd_t *md, turi_t *t, ttm_t *ttm)
{
    bool CalcSource = fabs(md->time - ttm->laser_t0) < ttm->laser_tdiff ? true : false;
    fmd_real_t dt = md->timestep / ttm->timestep_ratio;

    /* time loop */
    for (int j=0; j < ttm->timestep_ratio; j++)
    {
        fmd_real_t source_spcind;

        if (CalcSource)
            source_spcind = ttm->laser_factor_constant *
                            exp(ttm->laser_m_2sig2_inv * sqrr(md->time + j * dt - ttm->laser_t0));

        /* spatial loop */
        for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
        {
            if (ttm->Te_1d[i] < 0) continue;  /* nothing to be done with a deactivated cell */

            fmd_real_t source;

            if (CalcSource)
            {
                fmd_real_t z = (i - t->itc_glob_to_loc[2] + 0.5) * t->tcellh[2];
                source = source_spcind * exp(ttm->laser_m_absdepth_inv * (z - ttm->frontsurf));
            }
            else
                source = 0.0;

            int im1 = (i == t->itc_start_owned[2]) ? 0 : i-1;
            int ilo = (ttm->Te_1d[im1] < 0) ? i : im1;
            int ihi = (ttm->Te_1d[i+1] < 0) ? i : i+1;

            fmd_real_t KeL = 0.5 * (ttm->Ke_1d[i] + ttm->Ke_1d[ilo]);
            fmd_real_t KeR = 0.5 * (ttm->Ke_1d[i] + ttm->Ke_1d[ihi]);

            ttm->Te2_1d[i] = ttm->Te_1d[i] + dt/ttm->Ce_1d[i] * (
              (KeR * (ttm->Te_1d[ihi]-ttm->Te_1d[i]) - KeL * (ttm->Te_1d[i]-ttm->Te_1d[ilo])) / ttm->dz2
              - ttm->G_1d[i] * (ttm->Te_1d[i] - ttm->Ti_1d[i]) + source );

            ttm->xi_1d[i] += ttm->Te2_1d[i];
        } /* end of spatial loop */

        fmd_real_t *tempo = ttm->Te2_1d;
        ttm->Te2_1d = ttm->Te_1d;
        ttm->Te_1d = tempo;

        if (j < ttm->timestep_ratio-1)
        {
            type2_update_Ce_Ke_G(t, ttm);

            if (t->has_upper_lower_owner_procs[2]) _fmd_turi_update_ghosts_1d(md, t, 2, ttm->tgp);
        }
    } /* end of time loop */

    if (ttm->Te2_1d != ((fmd_real_t ***)ttm->Te_aux.data)[0][0])
    {
        fmd_array3s_t tempo = ttm->Te_aux;
        ttm->Te_aux = t->fields[ttm->iTe].data;
        t->fields[ttm->iTe].data = tempo;
    }

    for (int i = t->itc_start_owned[2]; i < t->itc_stop[2]; i++)
    {
        if (ttm->num_1d[i] >= ttm->min_atoms_num)
            ttm->xi_1d[i] = ttm->G_1d[i] * t->tcell_volume *
                            (ttm->xi_1d[i] / ttm->timestep_ratio - ttm->Ti_1d[i]) /
                            (3 * K_BOLTZMANN * ttm->num_1d[i] * ttm->Ti_1d[i]);
    }
}

static void init_common_types_1_2(fmd_t *md, ttm_t *ttm, turi_t *t)
{
    int inum = _fmd_field_add(t, FMD_FIELD_NUMBER, md->timestep, false);
    int ivcm = _fmd_field_add(t, FMD_FIELD_VCM, md->timestep, true);
    int iTi = _fmd_field_add(t, FMD_FIELD_TEMPERATURE, md->timestep, false);
    ttm->iTe = _fmd_field_add(t, FMD_FIELD_TTM_TE, md->timestep, false);
    ttm->ixi = _fmd_field_add(t, FMD_FIELD_TTM_XI, md->timestep, true);

    _fmd_array_3d_create(t->tdims_local, sizeof(fmd_real_t), DATATYPE_REAL, &ttm->Te_aux);

    /* set default value for "cell activation fraction" */
    ttm->CellActivFrac = 0.1;

    if (t->tdims_global[0] == 1 && t->tdims_global[1] == 1)
    {
        ttm->dim = 1;

        ttm->num_1d = ((unsigned ***)t->fields[inum].data.data)[0][0];
        ttm->vcm_1d = ((fmd_rtuple_t ***)t->fields[ivcm].data.data)[0][0];
        ttm->Ti_1d = ((fmd_real_t ***)t->fields[iTi].data.data)[0][0];
        ttm->Te_1d = ((fmd_real_t ***)t->fields[ttm->iTe].data.data)[0][0];
        ttm->Te2_1d = ((fmd_real_t ***)ttm->Te_aux.data)[0][0];
        ttm->xi_1d = ((fmd_real_t ***)t->fields[ttm->ixi].data.data)[0][0];

        ttm->num_1d[0] = ttm->num_1d[t->itc_stop[2]] = 0;

        ttm->dz2 = t->tcellh[2] * t->tcellh[2];

        if (t->ownerscomm.owned_tcells_num == 0)
            ttm->tgp = NULL;
        else
            ttm->tgp = (tghost_pack_t *)m_alloc(sizeof(tghost_pack_t));
    }
    else
    {
        ttm->dim = 3;

        ttm->num = (unsigned ***)t->fields[inum].data.data;
        ttm->vcm = (fmd_rtuple_t ***)t->fields[ivcm].data.data;
        ttm->Ti = (fmd_real_t ***)t->fields[iTi].data.data;
        ttm->Te = (fmd_real_t ***)t->fields[ttm->iTe].data.data;
        ttm->Te2 = (fmd_real_t ***)ttm->Te_aux.data;
        ttm->xi = (fmd_real_t ***)t->fields[ttm->ixi].data.data;
    }
}

static void ttm_init_type1(fmd_t *md, ttm_t *ttm, turi_t *t)
{
    init_common_types_1_2(md, ttm, t);

    if (ttm->dim == 1)
    {
        ttm->preupdate_xe_te = ttm_type1_presolve_1d;
        ttm->update_xe_te = ttm_type1_solve_1d;

        if (t->ownerscomm.owned_tcells_num > 0)
        {
            ttm->tgp->pack._1D = type1_1d_pack;
            ttm->tgp->unpack._1D = type1_1d_unpack;
            ttm->tgp->bufsize = type1_1d_calc_tghost_buffsize(md);
            ttm->tgp->sendbuf = m_alloc(ttm->tgp->bufsize);
            ttm->tgp->recvbuf = m_alloc(ttm->tgp->bufsize);
        }
    }
    else
    {
        ttm->preupdate_xe_te = NULL;
        ttm->update_xe_te = NULL; /* no 3D preupdater and updater is written */
    }
}

static void ttm_init_type2(fmd_t *md, ttm_t *ttm, turi_t *t)
{
    init_common_types_1_2(md, ttm, t);

    _fmd_array_3d_create(t->tdims_local, sizeof(fmd_real_t), DATATYPE_REAL, &ttm->Kel);
    _fmd_array_3d_create(t->tdims_local, sizeof(fmd_real_t), DATATYPE_REAL, &ttm->Cel);
    _fmd_array_3d_create(t->tdims_local, sizeof(fmd_real_t), DATATYPE_REAL, &ttm->Geis);

    ttm->T_C = ttm->Ct = ttm->C_DD = NULL;
    ttm->T_G = ttm->Gt = ttm->G_DD = NULL;

    if (ttm->dim == 1)
    {
        ttm->Ke_1d = ((fmd_real_t ***)ttm->Kel.data)[0][0];
        ttm->Ce_1d = ((fmd_real_t ***)ttm->Cel.data)[0][0];
        ttm->G_1d = ((fmd_real_t ***)ttm->Geis.data)[0][0];

        ttm->preupdate_xe_te = ttm_type2_presolve_1d;
        ttm->update_xe_te = ttm_type2_solve_1d;

        if (t->ownerscomm.owned_tcells_num > 0)
        {
            ttm->tgp->pack._1D = type2_1d_pack;
            ttm->tgp->unpack._1D = type2_1d_unpack;
            ttm->tgp->bufsize = type2_1d_calc_tghost_buffsize(md);
            ttm->tgp->sendbuf = m_alloc(ttm->tgp->bufsize);
            ttm->tgp->recvbuf = m_alloc(ttm->tgp->bufsize);
        }
    }
    else
    {
        ttm->Ke = (fmd_real_t ***)ttm->Kel.data;
        ttm->Ce = (fmd_real_t ***)ttm->Cel.data;
        ttm->Gei = (fmd_real_t ***)ttm->Geis.data;

        ttm->preupdate_xe_te = NULL;
        ttm->update_xe_te = NULL; /* no 3D preupdater and updater is written */
    }
}

ttm_t *_fmd_ttm_construct(fmd_t *md, turi_t *t)
{
    ttm_t *ttm = (ttm_t *)m_alloc(sizeof(ttm_t));

    switch (t->cat)
    {
        case FMD_TURI_TTM_TYPE1:
            ttm_init_type1(md, ttm, t);
            break;

        case FMD_TURI_TTM_TYPE2:
            ttm_init_type2(md, ttm, t);
            break;

        default:
            assert(0);
    }

    return ttm;
}

void _fmd_ttm_destruct(turi_t *t)
{
    _fmd_array_3d_free(&t->ttm->Te_aux);

    if (t->cat == FMD_TURI_TTM_TYPE2)
    {
        _fmd_array_3d_free(&t->ttm->Kel);
        _fmd_array_3d_free(&t->ttm->Cel);
        _fmd_array_3d_free(&t->ttm->Geis);
        free(t->ttm->T_C);
        free(t->ttm->C_DD);
        free(t->ttm->Ct);
        free(t->ttm->T_G);
        free(t->ttm->G_DD);
        free(t->ttm->Gt);
    }

    if (t->ttm->tgp != NULL)
    {
        free(t->ttm->tgp->recvbuf);
        free(t->ttm->tgp->sendbuf);
        free(t->ttm->tgp);
    }

    free(t->ttm);
    t->ttm = NULL;
}

void _fmd_ttm_setHeatCapacity_file(fmd_t *md, fmd_handle_t turi, fmd_string_t path)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    if (ttm->T_C != NULL)
    {
        free(ttm->T_C);
        free(ttm->Ct);
        free(ttm->C_DD);
        ttm->T_C = ttm->Ct = NULL;
    }

    FILE *fp = f_open(path, "r");

    double T, C;
    const int incr = 100;
    int i = 0, cap = 0;

    while (fscanf(fp, "%lf%lf", &T, &C) == 2)
    {
        if (i == cap)
        {
            cap += incr;
            ttm->T_C = re_alloc(ttm->T_C, cap * sizeof(fmd_real_t));
            ttm->Ct = re_alloc(ttm->Ct, cap * sizeof(fmd_real_t));
        }

        ttm->T_C[i] = T * 1e4;
        ttm->Ct[i] = C * 1e5 * JOULE_PER_METER3_KELVIN;
        i++;
    }

    fclose(fp);

    ttm->nC = i;

    if (ttm->nC < cap)
    {
        ttm->T_C = re_alloc(ttm->T_C, ttm->nC * sizeof(fmd_real_t));
        ttm->Ct = re_alloc(ttm->Ct, ttm->nC * sizeof(fmd_real_t));
    }

#ifdef USE_CSPLINE

    ttm->C_DD = m_alloc(ttm->nC * sizeof(fmd_real_t));
    spline_prepare_adv(ttm->T_C, ttm->Ct, ttm->nC, ttm->C_DD);

#endif

}

void _fmd_ttm_setHeatCapacity_linear(fmd_t *md, fmd_handle_t turi, fmd_ttm_heat_capacity_linear_t c)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    ttm->C_gamma = c.gamma * JOULE_PER_METER3_KELVIN2;
}

void _fmd_ttm_setHeatConductivity_zhigilei(fmd_t *md, fmd_handle_t turi, fmd_ttm_heat_conductivity_zhigilei_t k)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    ttm->Kzh.A = k.A / SECOND;
    ttm->Kzh.B = k.B / SECOND;
    ttm->Kzh.v2 = sqrr(k.v * METER_PER_SECOND);
}

void _fmd_ttm_setHeatConductivity_constant1(fmd_t *md, fmd_handle_t turi, fmd_real_t k)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    ttm->K = k * WATT_PER_METER_KELVIN;
}

void _fmd_ttm_setHeatConductivity_constant2(fmd_t *md, fmd_handle_t turi, fmd_ttm_heat_conductivity_constant_t k)
{
    _fmd_ttm_setHeatConductivity_constant1(md, turi, k.value);
}

void _fmd_ttm_setCouplingFactor_file(fmd_t *md, fmd_handle_t turi, fmd_string_t path)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    if (ttm->T_G != NULL)
    {
        free(ttm->T_G);
        free(ttm->Gt);
        free(ttm->G_DD);
        ttm->T_G = ttm->Gt = NULL;
    }

    FILE *fp = f_open(path, "r");

    double T, G;
    const int incr = 100;
    int i = 0, cap = 0;

    while (fscanf(fp, "%lf%lf", &T, &G) == 2)
    {
        if (i == cap)
        {
            cap += incr;
            ttm->T_G = re_alloc(ttm->T_G, cap * sizeof(fmd_real_t));
            ttm->Gt = re_alloc(ttm->Gt, cap * sizeof(fmd_real_t));
        }

        ttm->T_G[i] = T * 1e4;
        ttm->Gt[i] = G * 1e17 * WATT_PER_METER3_KELVIN;
        i++;
    }

    fclose(fp);

    ttm->nG = i;

    if (ttm->nG < cap)
    {
        ttm->T_G = re_alloc(ttm->T_G, ttm->nG * sizeof(fmd_real_t));
        ttm->Gt = re_alloc(ttm->Gt, ttm->nG * sizeof(fmd_real_t));
    }

#ifdef USE_CSPLINE

    ttm->G_DD = m_alloc(ttm->nG * sizeof(fmd_real_t));
    spline_prepare_adv(ttm->T_G, ttm->Gt, ttm->nG, ttm->G_DD);

#endif

}

void _fmd_ttm_setCouplingFactor_constant1(fmd_t *md, fmd_handle_t turi, fmd_real_t g)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    ttm->G = g * WATT_PER_METER3_KELVIN;
}

void _fmd_ttm_setCouplingFactor_constant2(fmd_t *md, fmd_handle_t turi, fmd_ttm_coupling_factor_constant_t g)
{
    _fmd_ttm_setCouplingFactor_constant1(md, turi, g.value);
}

void fmd_ttm_setElectronTemperature(fmd_t *md, fmd_handle_t turi, fmd_real_t Te)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1 || t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    fmd_ituple_t itc;

    LOOP3D(itc, t->itc_start, t->itc_stop)
        ARRAY_ELEMENT((fmd_real_t ***)t->fields[ttm->iTe].data.data, itc) = Te;
}

void fmd_ttm_setTimestepRatio(fmd_t *md, fmd_handle_t turi, int ratio)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1 || t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    ttm->timestep_ratio = ratio;
}

void fmd_ttm_setCellActivationFraction(fmd_t *md, fmd_handle_t turi, fmd_real_t value)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1 || t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    assert(value >= 0 && value <= 1); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    ttm->CellActivFrac = value;
}

void _fmd_ttm_setLaserSource_gaussian(fmd_t *md, fmd_handle_t turi, fmd_ttm_laser_gaussian_t laser)
{
    turi_t *t = &md->turies[turi];

    assert(t->cat == FMD_TURI_TTM_TYPE1 || t->cat == FMD_TURI_TTM_TYPE2); /* TO-DO: handle error */

    ttm_t *ttm = t->ttm;

    fmd_real_t Lp = laser.AbsorptionDepth * METER;
    fmd_real_t R = laser.reflectance;
    fmd_real_t I0 = sqrt(4*log(2)/M_PI) * laser.fluence / laser.duration * WATT_PER_METER2;
    fmd_real_t sigma = laser.duration / sqrt(8*log(2)) * SECOND;

    ttm->laser_t0 = laser.t0 * SECOND;
    ttm->laser_tdiff = 4*sqrt(2) * sigma;
    ttm->laser_m_2sig2_inv = -0.5 / sqrr(sigma);
    ttm->laser_m_absdepth_inv = -1. / Lp;
    ttm->laser_factor_constant = I0 * (1 - R) / Lp;
}
