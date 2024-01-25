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
#define EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum(c1, i1, atomkind1, rho_host, atomkinds,  \
                                                  Femb_sum)                                \
    do                                                                                     \
    {                                                                                      \
        eam_element_t *el = atomkinds[atomkind1].eam_element;                              \
        eam_t *eam = el->eam;                                                              \
        fmd_real_t h = eam->drho;                                                          \
        int irho = (int)(rho_host / h);                                                    \
        assert(irho < eam->Nrho - 1);                                                      \
        int irho_h = irho + 1;                                                             \
        fmd_real_t *F = el->F;                                                             \
        fmd_real_t *F_DD = el->F_DD;                                                       \
        fmd_real_t a = irho_h - rho_host/h;                                                \
        fmd_real_t b = 1-a;                                                                \
        c1->FembPrime[i1] = SPLINE_DERIV(a,b,F,irho,irho_h,F_DD,h);                        \
        Femb_sum += SPLINE_VAL(a,b,F,irho,irho_h,F_DD,h);                                  \
    } while (0)
#else
#define EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum(c1, i1, atomkind1, rho_host, atomkinds,  \
                                                  Femb_sum)                                \
    do                                                                                     \
    {                                                                                      \
        eam_element_t *el = atomkinds[atomkind1].eam_element;                              \
        eam_t *eam = el->eam;                                                              \
        fmd_real_t h = eam->drho;                                                          \
        int irho = (int)(rho_host / h);                                                    \
        assert(irho < eam->Nrho - 1);                                                      \
        int irho_h = irho + 1;                                                             \
        fmd_real_t *F = el->F;                                                             \
        c1->FembPrime[i1] = (F[irho_h] - F[irho]) / h;                                     \
        Femb_sum += F[irho] + (rho_host - irho * h) * c1->FembPrime[i1];                   \
    } while (0)
#endif

#ifdef USE_CSPLINE
#define EAM_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, kc, c1, i1, atomkind2,          \
                                            c2, i2, PotEn, pottable)                       \
    do                                                                                     \
    {                                                                                      \
        fmd_real_t r2;                                                                     \
        fmd_rtuple_t rv;                                                                   \
                                                                                           \
        fmd_real_t *x2 = &POS(c2, i2, 0);                                                  \
                                                                                           \
        COMPUTE_rv_AND_r2(x1, x2, kc, rv, r2);                                             \
                                                                                           \
        eam_t *eam = (eam_t *)pottable[atomkind1][atomkind2].data;                         \
                                                                                           \
        if (r2 < eam->cutoff_sqr)                                                          \
        {                                                                                  \
            fmd_real_t h = eam->dr2;                                                       \
            int ir2 = (int)(r2 / h);                                                       \
            int ir2_h = ir2 + 1;                                                           \
            unsigned iloc = pottable[atomkind1][atomkind2].iloc;                           \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;                           \
                                                                                           \
            fmd_real_t *rho_i = eam->elements[iloc].rho;                                   \
            fmd_real_t *rho_j;                                                             \
            fmd_real_t *rho_iDD = eam->elements[iloc].rhoDD;                               \
            fmd_real_t *rho_jDD;                                                           \
            fmd_real_t *phi = eam->elements[iloc].phi[jloc];                               \
            fmd_real_t *phiDD = eam->elements[iloc].phiDD[jloc];                           \
            fmd_real_t a = ir2_h - r2/h;                                                   \
            fmd_real_t b = 1-a;                                                            \
            fmd_real_t phi_deriv = SPLINE_DERIV(a,b,phi,ir2,ir2_h,phiDD,h);                \
            fmd_real_t rho_ip = SPLINE_DERIV(a,b,rho_i,ir2,ir2_h,rho_iDD,h);               \
            fmd_real_t rho_jp;                                                             \
                                                                                           \
            if (jloc == iloc)                                                              \
                rho_jp = rho_ip;                                                           \
            else                                                                           \
            {                                                                              \
                rho_j = eam->elements[jloc].rho;                                           \
                rho_jDD = eam->elements[jloc].rhoDD;                                       \
                rho_jp = SPLINE_DERIV(a,b,rho_j,ir2,ir2_h,rho_jDD,h);                      \
            }                                                                              \
                                                                                           \
            fmd_real_t mag = 2 * (c1->FembPrime[i1] * rho_jp +                             \
                                  c2->FembPrime[i2] * rho_ip + phi_deriv);                 \
            PotEn += SPLINE_VAL(a,b,phi,ir2,ir2_h,phiDD,h);                                \
                                                                                           \
            for (int d=0; d<DIM; d++)                                                      \
                FRC(c1, i1, d) -= mag * rv[d];                                             \
        }                                                                                  \
    } while (0)
#else
#define EAM_PAIR_UPDATE_FORCE_AND_POTENERGY(x1, atomkind1, kc, c1, i1, atomkind2,          \
                                            c2, i2, PotEn, pottable)                       \
    do                                                                                     \
    {                                                                                      \
        fmd_real_t r2;                                                                     \
        fmd_rtuple_t rv;                                                                   \
                                                                                           \
        fmd_real_t *x2 = &POS(c2, i2, 0);                                                  \
                                                                                           \
        COMPUTE_rv_AND_r2(x1, x2, kc, rv, r2);                                             \
                                                                                           \
        eam_t *eam = (eam_t *)pottable[atomkind1][atomkind2].data;                         \
                                                                                           \
        if (r2 < eam->cutoff_sqr)                                                          \
        {                                                                                  \
            fmd_real_t h = eam->dr2;                                                       \
            int ir2 = (int)(r2 / h);                                                       \
            int ir2_h = ir2 + 1;                                                           \
            unsigned iloc = pottable[atomkind1][atomkind2].iloc;                           \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;                           \
                                                                                           \
            fmd_real_t *rho_i = eam->elements[iloc].rho;                                   \
            fmd_real_t *phi = eam->elements[iloc].phi[jloc];                               \
            fmd_real_t *rho_j = eam->elements[jloc].rho;                                   \
            fmd_real_t mag = 2 * (c1->FembPrime[i1] * (rho_j[ir2_h] - rho_j[ir2]) +        \
                                  c2->FembPrime[i2] * (rho_i[ir2_h] - rho_i[ir2]) +        \
                                  (phi[ir2_h] - phi[ir2])) / h;                            \
            PotEn += phi[ir2] + (r2/h - ir2) * (phi[ir2_h] - phi[ir2]);                    \
                                                                                           \
            for (int d=0; d<DIM; d++)                                                      \
                FRC(c1, i1, d) -= mag * rv[d];                                             \
        }                                                                                  \
    } while (0)
#endif

#ifdef USE_CSPLINE
#define EAM_PAIR_UPDATE_rho_host(x1, atomkind1, kc, c1, i1, atomkind2, c2, i2, rho_host,   \
                                 pottable)                                                 \
    do                                                                                     \
    {                                                                                      \
        fmd_real_t r2;                                                                     \
        fmd_rtuple_t rv;                                                                   \
                                                                                           \
        fmd_real_t *x2 = &POS(c2, i2, 0);                                                  \
                                                                                           \
        COMPUTE_rv_AND_r2(x1, x2, kc, rv, r2);                                             \
                                                                                           \
        eam_t *eam = (eam_t *)pottable[atomkind1][atomkind2].data;                         \
                                                                                           \
        if (r2 < eam->cutoff_sqr)                                                          \
        {                                                                                  \
            fmd_real_t h = eam->dr2;                                                       \
            int ir2 = (int)(r2 / h);                                                       \
            int ir2_h = ir2 + 1;                                                           \
            fmd_real_t a = ir2_h - r2/h;                                                   \
            fmd_real_t b=1-a;                                                              \
                                                                                           \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;                           \
            fmd_real_t *rho = eam->elements[jloc].rho;                                     \
            fmd_real_t *rhoDD = eam->elements[jloc].rhoDD;                                 \
            rho_host += SPLINE_VAL(a,b,rho,ir2,ir2_h,rhoDD,h);                             \
        }                                                                                  \
    } while (0)
#else
#define EAM_PAIR_UPDATE_rho_host(x1, atomkind1, kc, c1, i1, atomkind2, c2, i2, rho_host,   \
                                 pottable)                                                 \
    do                                                                                     \
    {                                                                                      \
        fmd_real_t r2;                                                                     \
        fmd_rtuple_t rv;                                                                   \
                                                                                           \
        fmd_real_t *x2 = &POS(c2, i2, 0);                                                  \
                                                                                           \
        COMPUTE_rv_AND_r2(x1, x2, kc, rv, r2);                                             \
                                                                                           \
        eam_t *eam = (eam_t *)pottable[atomkind1][atomkind2].data;                         \
                                                                                           \
        if (r2 < eam->cutoff_sqr)                                                          \
        {                                                                                  \
            fmd_real_t h = eam->dr2;                                                       \
            int ir2 = (int)(r2 / h);                                                       \
            int ir2_h = ir2 + 1;                                                           \
            fmd_real_t a = ir2_h - r2/h;                                                   \
            fmd_real_t b = 1 - a;                                                          \
                                                                                           \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;                           \
            fmd_real_t *rho = eam->elements[jloc].rho;                                     \
            rho_host += rho[ir2]*a + rho[ir2_h]*b;                                         \
        }                                                                                  \
    } while (0)
#endif

typedef struct _eam eam_t;

typedef struct _eam_element
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

struct _eam
{
    eam_element_t *elements;
    fmd_real_t drho, dr, dr2, cutoff_sqr;
    int ElementsNo;
    int Nrho, Nr, Nr2;
};

void _fmd_computeEAM_pass0(fmd_t *md, fmd_real_t FembSum);
void _fmd_computeEAM_pass1(fmd_t *md, fmd_real_t *FembSum_p);
unsigned _fmd_pot_eam_find_iloc(fmd_t *md, eam_t *eam, unsigned atomkind);
void _fmd_pot_eam_free(eam_t *eam);

#endif /* EAM_H */
