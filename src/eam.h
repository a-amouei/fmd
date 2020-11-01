/*
  eam.h: This file is part of Free Molecular Dynamics

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

#ifndef EAM_H
#define EAM_H

#include "config.h"
#include "forces.h"
#include "types.h"
#include "cspline.h"

#ifdef USE_CSPLINE
#define EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum                            \
    {                                                                        \
        eam_element_t *el = atomkinds[atomkind1].eam_element;                \
        eam = el->eam;                                                       \
        h = eam->drho;                                                       \
        irho = (int)(rho_host / h);                                          \
        assert(irho < eam->Nrho - 1);                                        \
        irho_h = irho + 1;                                                   \
        F = el->F;                                                           \
        F_DD = el->F_DD;                                                     \
        a = irho_h - rho_host/h;                                             \
        b = 1-a;                                                             \
        item1_p->FembPrime = SPLINE_DERIV(a,b,F,irho,irho_h,F_DD,h);         \
        Femb_sum += SPLINE_VAL(a,b,F,irho,irho_h,F_DD,h);                    \
    }
#else
#define EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum                            \
    {                                                                        \
        eam_element_t *el = atomkinds[atomkind1].eam_element;                \
        eam = el->eam;                                                       \
        h = eam->drho;                                                       \
        irho = (int)(rho_host / h);                                          \
        assert(irho < eam->Nrho - 1);                                        \
        irho_h = irho + 1;                                                   \
        F = el->F;                                                           \
        item1_p->FembPrime = (F[irho_h] - F[irho]) / h;                      \
        Femb_sum += F[irho] + (rho_host - irho * h) * item1_p->FembPrime;    \
    }
#endif

#ifdef USE_CSPLINE
#define EAM_PAIR_UPDATE_FORCE_AND_POTENERGY                                     \
    {                                                                           \
        COMPUTE_rv_AND_r2;                                                      \
                                                                                \
        eam = (eam_t *)pottable[atomkind1][atomkind2].data;                     \
                                                                                \
        if (r2 < eam->cutoff_sqr)                                               \
        {                                                                       \
            h = eam->dr2;                                                       \
            ir2 = (int)(r2 / h);                                                \
            ir2_h = ir2 + 1;                                                    \
            unsigned iloc = pottable[atomkind1][atomkind2].iloc;                \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;                \
                                                                                \
            rho_i = eam->elements[iloc].rho;                                    \
            rho_iDD = eam->elements[iloc].rhoDD;                                \
            phi = eam->elements[iloc].phi[jloc];                                \
            phiDD = eam->elements[iloc].phiDD[jloc];                            \
            a = ir2_h - r2/h;                                                   \
            b = 1-a;                                                            \
            phi_deriv = SPLINE_DERIV(a,b,phi,ir2,ir2_h,phiDD,h);                \
            rho_ip = SPLINE_DERIV(a,b,rho_i,ir2,ir2_h,rho_iDD,h);               \
            if (jloc == iloc)                                                   \
                rho_jp = rho_ip;                                                \
            else                                                                \
            {                                                                   \
                rho_j = eam->elements[jloc].rho;                                \
                rho_jDD = eam->elements[jloc].rhoDD;                            \
                rho_jp = SPLINE_DERIV(a,b,rho_j,ir2,ir2_h,rho_jDD,h);           \
            }                                                                   \
                                                                                \
            mag = 2 * (item1_p->FembPrime * rho_jp +                            \
                  item2_p->FembPrime * rho_ip + phi_deriv);                     \
            potEnergy += SPLINE_VAL(a,b,phi,ir2,ir2_h,phiDD,h);                 \
                                                                                \
            for (d=0; d<3; d++)                                                 \
                item1_p->F[d] -= mag * rv[d];                                   \
        }                                                                       \
    }
#else
#define EAM_PAIR_UPDATE_FORCE_AND_POTENERGY                                     \
    {                                                                           \
        COMPUTE_rv_AND_r2;                                                      \
                                                                                \
        eam = (eam_t *)pottable[atomkind1][atomkind2].data;                     \
                                                                                \
        if (r2 < eam->cutoff_sqr)                                               \
        {                                                                       \
            h = eam->dr2;                                                       \
            ir2 = (int)(r2 / h);                                                \
            ir2_h = ir2 + 1;                                                    \
            unsigned iloc = pottable[atomkind1][atomkind2].iloc;                \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;                \
                                                                                \
            rho_i = eam->elements[iloc].rho;                                    \
            phi = eam->elements[iloc].phi[jloc];                                \
            rho_j = eam->elements[jloc].rho;                                    \
            mag = 2 * (item1_p->FembPrime * (rho_j[ir2_h] - rho_j[ir2]) +       \
                       item2_p->FembPrime * (rho_i[ir2_h] - rho_i[ir2]) +       \
                                                (phi[ir2_h] - phi[ir2])) / h;   \
            potEnergy += phi[ir2] + (r2/h - ir2) * (phi[ir2_h] - phi[ir2]);     \
                                                                                \
            for (d=0; d<3; d++)                                                 \
                item1_p->F[d] -= mag * rv[d];                                   \
        }                                                                       \
    }
#endif

#ifdef USE_CSPLINE
#define EAM_PAIR_UPDATE_rho_host                                       \
    {                                                                  \
        COMPUTE_rv_AND_r2;                                             \
                                                                       \
        eam = (eam_t *)pottable[atomkind1][atomkind2].data;            \
                                                                       \
        if (r2 < eam->cutoff_sqr)                                      \
        {                                                              \
            h = eam->dr2;                                              \
            ir2 = (int)(r2 / h);                                       \
            ir2_h = ir2 + 1;                                           \
            a = ir2_h - r2/h;                                          \
            b=1-a;                                                     \
                                                                       \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;       \
            rho = eam->elements[jloc].rho;                             \
            rhoDD = eam->elements[jloc].rhoDD;                         \
            rho_host += SPLINE_VAL(a,b,rho,ir2,ir2_h,rhoDD,h);         \
        }                                                              \
    }
#else
#define EAM_PAIR_UPDATE_rho_host                                       \
    {                                                                  \
        COMPUTE_rv_AND_r2;                                             \
                                                                       \
        eam = (eam_t *)pottable[atomkind1][atomkind2].data;            \
                                                                       \
        if (r2 < eam->cutoff_sqr)                                      \
        {                                                              \
            h = eam->dr2;                                              \
            ir2 = (int)(r2 / h);                                       \
            ir2_h = ir2 + 1;                                           \
            a = ir2_h - r2/h;                                          \
            b=1-a;                                                     \
                                                                       \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;       \
            rho = eam->elements[jloc].rho;                             \
            rho_host += rho[ir2]*a + rho[ir2_h]*b;                     \
        }                                                              \
    }
#endif

typedef struct eam_t eam_t;

typedef struct eam_element_t
{
    fmd_real_t mass;
    fmd_real_t latticeParameter;
    fmd_real_t *F;
    fmd_real_t *F_DD;
    fmd_real_t *rho;
    fmd_real_t *rhoDD;
    fmd_real_t **phi;
    fmd_real_t **phiDD;
    fmd_string_t name;
    eam_t *eam;
} eam_element_t;

struct eam_t
{
    eam_element_t *elements;
    fmd_real_t drho, dr, dr2, cutoff_sqr;
    int elementsNo;
    int Nrho, Nr, Nr2;
};

void fmd_computeEAM_pass0(fmd_t *md, fmd_real_t FembSum);
void fmd_computeEAM_pass1(fmd_t *md, fmd_real_t *FembSum_p);
unsigned fmd_pot_eam_find_iloc(fmd_t *md, eam_t *eam, unsigned atomkind);
void fmd_pot_eam_free(eam_t *eam);

#endif /* EAM_H */
