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
#include "molecule.h"

static void compute_hybrid_pass1(fmd_t *md, fmd_real_t *FembSum_p)
{
    fmd_ituple_t jc, kc;
    int d, ir2, irho, ir2_h, irho_h;
    fmd_real_t r2;
    fmd_rtuple_t rv;
    fmd_real_t *rho, *rhoDD, *F, *F_DD;
    fmd_real_t a, b, h;
    int ic0, ic1, ic2;
    fmd_real_t Femb_sum=0.0;
    potpair_t **pottable = md->potsys.pottable;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    /* iterate over all cells */
    #pragma omp parallel for private(ic0,ic1,ic2,kc,jc,d,rv,r2,h,ir2,ir2_h,a,b,rho, \
      rhoDD,F,F_DD,irho,irho_h) shared(md,pottable,atomkinds) default(none) collapse(3) reduction(+:Femb_sum) \
      schedule(static,1)
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
    for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
    for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
    {
        /* iterate over all particles in cell ic */
        cell_t *c1;
        int pind1;

        for (c1 = &md->SubDomain.grid[ic0][ic1][ic2], pind1=0; pind1 < c1->parts_num; pind1++)
        {
            particle_t *p1 = &c1->parts[pind1];

            if (!(md->ActiveGroup == -1 || p1->core.GroupID == md->ActiveGroup))
                continue;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = p1->core.atomkind;

            if (md->potsys.atomkinds[atomkind1].eam_element == NULL)
                continue;

            fmd_real_t rho_host = 0.0;

            /* iterate over neighbor cells of cell ic */
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
                        cell_t *c2;
                        int pind2;

                        for (c2 = &md->SubDomain.grid[jc[0]][jc[1]][jc[2]], pind2=0; pind2 < c2->parts_num; pind2++)
                        {
                            particle_t *p2 = &c2->parts[pind2];

                            if (!(md->ActiveGroup == -1 || p2->core.GroupID == md->ActiveGroup))
                                continue;

                            if (p1 != p2)
                            {
                                atomkind2 = p2->core.atomkind;

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

static void compute_hybrid_pass0(fmd_t *md, fmd_real_t FembSum)
{
    fmd_ituple_t jc, kc;
    int d, ir2, ir2_h;
    fmd_real_t r2;
    fmd_rtuple_t rv;
    fmd_real_t *rho_i, *rho_j, *phi;
    fmd_real_t *rho_iDD, *rho_jDD, *phiDD;
    fmd_real_t rho_ip, rho_jp;
    fmd_real_t mag;
    fmd_real_t phi_deriv;
    fmd_real_t a, b, h;
    int ic0, ic1, ic2;
    potpair_t **pottable = md->potsys.pottable;
#ifdef USE_TTM
    fmd_real_t mass;
    int ttm_index;
    fmd_real_t dx;
    fmd_real_t pxx = 0.0;
#endif
    fmd_real_t potEnergy = 0.0;

    /* iterate over all cells */
#ifdef USE_TTM
    #pragma omp parallel for private(ic0,ic1,ic2,ttm_index,d,element_i,rho_i,rho_iDD,kc,jc,rv,r2,h,ir2, \
      ir2_h,element_j,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag,mass,dx) \
      shared(md,ttm_lattice_aux,ttm_useSuction,ttm_suctionWidth,ttm_suctionIntensity,ttm_pxx_compute, \
      ttm_pxx_pos) default(none) collapse(3) reduction(+:potEnergy,pxx) schedule(static,1)
#else
    #pragma omp parallel for private(ic0,ic1,ic2,d,rho_i,rho_iDD,kc,jc,rv,r2,h,ir2, \
      ir2_h,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag) \
      shared(md,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
#endif
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
    for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
    for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
    {
#ifdef USE_TTM
        ttm_index = ic0 - md->SubDomain.ic_start[0] + 1;
#endif

        /* iterate over all particles in cell ic */

        cell_t *c1;
        int pind1;

        for (c1 = &md->SubDomain.grid[ic0][ic1][ic2], pind1=0; pind1 < c1->parts_num; pind1++)
        {
            particle_t *p1 = &c1->parts[pind1];

            if (!(md->ActiveGroup == -1 || p1->core.GroupID == md->ActiveGroup))
                continue;

            for (d=0; d<3; d++)
                p1->F[d] = 0.0;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = p1->core.atomkind;

            /* iterate over neighbor cells of cell ic */
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)

                        /* iterate over all items in cell jc */

                        cell_t *c2;
                        int pind2;

                        for (c2 = &md->SubDomain.grid[jc[0]][jc[1]][jc[2]], pind2=0; pind2 < c2->parts_num; pind2++)
                        {
                            particle_t *p2 = &c2->parts[pind2];

                            if (!(md->ActiveGroup == -1 || p2->core.GroupID == md->ActiveGroup))
                                continue;

                            if (p1 != p2)
                            {
                                atomkind2 = p2->core.atomkind;

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
    MPI_Allreduce(&potEnergy, &md->TotalPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
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
                fmd_real_t FembSum;
                fmd_computeEAM_pass1(md, &FembSum);
                fmd_ghostparticles_update_Femb(md);
                fmd_computeEAM_pass0(md, FembSum);
                break;
        }
    }
    else  // hybrid mode
    {
        fmd_real_t FembSum = 0.0;

        if (md->potsys.hybridpasses[1])
        {
            compute_hybrid_pass1(md, &FembSum);
            fmd_ghostparticles_update_Femb(md);
        }

        if (md->potsys.hybridpasses[0])
            compute_hybrid_pass0(md, FembSum);
    }

    if (md->TotalNoOfMolecules > 0) fmd_dync_computeBondForce(md);   // TO-DO

    fmd_ghostparticles_delete(md);
}
