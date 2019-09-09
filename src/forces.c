/*
  forces.c: This file is part of Free Molecular Dynamics

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

#include "forces.h"
#include "eam.h"
#include "lj.h"
#include "morse.h"
#include "base.h"
#include "md_ghost.h"
#include "list.h"

static void compute_hybrid_pass1(fmd_t *md, double *FembSum_p)
{
    int jc[3], kc[3];
    int d, ir2, irho, ir2_h, irho_h;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    double *rho, *rhoDD, *F, *F_DD;
    double a, b, h;
    int ic0, ic1, ic2;
    double Femb_sum=0.0;
    potpair_t **pottable = md->potsys.pottable;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,kc,jc,item2_p,d,rv,r2,h,ir2,ir2_h,a,b,rho, \
      rhoDD,F,F_DD,irho,irho_h) shared(md,pottable,atomkinds) default(none) collapse(3) reduction(+:Femb_sum) \
      schedule(static,1)
    for (ic0 = md->subDomain.ic_start[0]; ic0 < md->subDomain.ic_stop[0]; ic0++)
    for (ic1 = md->subDomain.ic_start[1]; ic1 < md->subDomain.ic_stop[1]; ic1++)
    for (ic2 = md->subDomain.ic_start[2]; ic2 < md->subDomain.ic_stop[2]; ic2++)
    {
        // iterate over all items in cell ic
        for (item1_p = md->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(md->activeGroup == -1 || item1_p->P.GroupID == md->activeGroup))
                continue;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = item1_p->P.atomkind;

            if (md->potsys.atomkinds[atomkind1].eam_element == NULL)
                continue;

            double rho_host = 0.0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = md->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(md->activeGroup == -1 || item2_p->P.GroupID == md->activeGroup))
                                continue;

                            if (item1_p != item2_p)
                            {
                                atomkind2 = item2_p->P.atomkind;
                                if (pottable[atomkind1][atomkind2].cat == POT_EAM_ALLOY)
                                    EAM_PAIR_UPDATE_rho_host;
                            }
                        }
                    }
                }
            }

            EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum;
        }
    }

    *FembSum_p=Femb_sum;
}

static void compute_hybrid_pass0(fmd_t *md, double FembSum)
{
    int jc[3], kc[3];
    int d, ir2, ir2_h;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    double *rho_i, *rho_j, *phi;
    double *rho_iDD, *rho_jDD, *phiDD;
    double rho_ip, rho_jp;
    double mag;
    double phi_deriv;
    double a, b, h;
    int ic0, ic1, ic2;
    potpair_t **pottable = md->potsys.pottable;
#ifdef USE_TTM
    double mass;
    int ttm_index;
    double dx;
    double pxx = 0.0;
#endif
    double potEnergy = 0.0;

    // iterate over all cells(lists)
#ifdef USE_TTM
    #pragma omp parallel for private(ic0,ic1,ic2,ttm_index,item1_p,d,element_i,rho_i,rho_iDD,kc,jc,item2_p,rv,r2,h,ir2, \
      ir2_h,element_j,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag,mass,dx) \
      shared(md,ttm_lattice_aux,ttm_useSuction,ttm_suctionWidth,ttm_suctionIntensity,ttm_pxx_compute, \
      ttm_pxx_pos) default(none) collapse(3) reduction(+:potEnergy,pxx) schedule(static,1)
#else
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,rho_i,rho_iDD,kc,jc,item2_p,rv,r2,h,ir2, \
      ir2_h,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag) \
      shared(md,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
#endif
    for (ic0 = md->subDomain.ic_start[0]; ic0 < md->subDomain.ic_stop[0]; ic0++)
    for (ic1 = md->subDomain.ic_start[1]; ic1 < md->subDomain.ic_stop[1]; ic1++)
    for (ic2 = md->subDomain.ic_start[2]; ic2 < md->subDomain.ic_stop[2]; ic2++)
    {
#ifdef USE_TTM
        ttm_index = ic0 - md->subDomain.ic_start[0] + 1;
#endif
        // iterate over all items in cell ic
        for (item1_p = md->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(md->activeGroup == -1 || item1_p->P.GroupID == md->activeGroup))
                continue;

            for (d=0; d<3; d++)
                item1_p->F[d] = 0.0;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = item1_p->P.atomkind;

            // iterate over neighbor cells of cell ic
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = md->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(md->activeGroup == -1 || item2_p->P.GroupID == md->activeGroup))
                                continue;

                            if (item1_p != item2_p)
                            {
                                atomkind2 = item2_p->P.atomkind;

                                switch (pottable[atomkind1][atomkind2].cat)
                                {
                                    case POT_EAM_ALLOY:
                                        EAM_PAIR_UPDATE_FORCE_AND_POTENERGY;
                                        break;

                                    case POT_LJ_6_12:
                                        LJ_PAIR_UPDATE_FORCE_AND_POTENERGY;
                                        break;

                                    case POT_MORSE:
                                        MORSE_PAIR_UPDATE_FORCE_AND_POTENERGY;
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    potEnergy = 0.5 * potEnergy + FembSum;
    MPI_Allreduce(&potEnergy, &md->totalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, md->MD_comm);
}

void fmd_dync_updateForces(fmd_t *md)
{
    if (md->potsys.potcats == NULL)  // just for one time
        fmd_pot_prepareForForceComp(md);

    fmd_ghostparticles_init(md);

    if (md->potsys.potcats_num == 1) // not hybrid mode
    {
        potcat_t potkind = *(potcat_t *)(md->potsys.potcats->data);

        switch (potkind)
        {
            case POT_LJ_6_12:
                fmd_computeLJ(md);
                break;

            case POT_MORSE:
                fmd_computeMorse(md);
                break;

            case POT_EAM_ALLOY:
                if (md->CompLocOrdParam) compLocOrdParam(md);
                double FembSum;
                fmd_computeEAM_pass1(md, &FembSum);
                fmd_ghostparticles_update_Femb(md);
                fmd_computeEAM_pass0(md, FembSum);
                break;
        }
    }
    else  // hybrid mode
    {
        double FembSum = 0.0;

        if (md->potsys.hybridpasses[1])
        {
            compute_hybrid_pass1(md, &FembSum);
            fmd_ghostparticles_update_Femb(md);
        }

        if (md->potsys.hybridpasses[0])
            compute_hybrid_pass0(md, FembSum);
    }

    fmd_ghostparticles_delete(md);
}
