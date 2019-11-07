/*
  eam.c: This file is part of Free Molecular Dynamics

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

#include "eam.h"
#include "base.h"
#include "list.h"

#define EAM_TTM_UPDATE_FORCE                                                    \
    {                                                                           \
        mass = md->EAM.elements[element_i].mass;                              \
        for (d=0; d<3; d++)                                                     \
            item1_p->F[d] += ttm_lattice_aux[ttm_index].xi *                    \
                mass * (item1_p->P.v[d] - ttm_lattice_aux[ttm_index].v_cm[d]);  \
        if (ttm_useSuction)                                                     \
            if (item1_p->P.x[0] < ttm_suctionWidth)                             \
                item1_p->F[0] -= mass * ttm_suctionIntensity;                   \
        if (ttm_pxx_compute)                                                    \
        {                                                                       \
            dx = item1_p->P.x[0] - ttm_pxx_pos;                                 \
            pxx += item1_p->F[0] * ((dx > 0) - (dx < 0));                       \
        }                                                                       \
    }

void fmd_computeEAM_pass0(fmd_t *md, double FembSum)
{
    int jc[3], kc[3];
    int d, ir2, ir2_h;
    ParticleListItem_t *item1_p, *item2_p;
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
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
    for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
    for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
    {
#ifdef USE_TTM
        ttm_index = ic0 - md->SubDomain.ic_start[0] + 1;
#endif
        // iterate over all items in cell ic
        for (item1_p = md->SubDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(md->ActiveGroup == -1 || item1_p->P.GroupID == md->ActiveGroup))
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
                        for (item2_p = md->SubDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(md->ActiveGroup == -1 || item2_p->P.GroupID == md->ActiveGroup))
                                continue;

                            if (item1_p != item2_p)
                            {
                                atomkind2 = item2_p->P.atomkind;
                                EAM_PAIR_UPDATE_FORCE_AND_POTENERGY;
                            }
                        }
                    }
                }
            }

#ifdef USE_TTM
            EAM_TTM_UPDATE_FORCE;
#endif
        }
    }

#ifdef USE_TTM
    ttm_pxx_local[1] += pxx;
#endif

    potEnergy = 0.5 * potEnergy + FembSum;
    MPI_Allreduce(&potEnergy, &md->TotalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, md->MD_comm);
}

void fmd_computeEAM_pass1(fmd_t *md, double *FembSum_p)
{
    int jc[3], kc[3];
    int d, ir2, irho, ir2_h, irho_h;
    ParticleListItem_t *item1_p, *item2_p;
    double r2, rv[3];
    double *rho, *rhoDD, *F, *F_DD;
    double a, b, h;
    int ic0, ic1, ic2;
    double Femb_sum=0;
    potpair_t **pottable = md->potsys.pottable;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,kc,jc,item2_p,d,rv,r2,h,ir2,ir2_h,a,b,rho, \
      rhoDD,F,F_DD,irho,irho_h) shared(md,pottable,atomkinds) default(none) collapse(3) reduction(+:Femb_sum) \
      schedule(static,1)
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
    for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
    for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
    {
        // iterate over all items in cell ic
        for (item1_p = md->SubDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(md->ActiveGroup == -1 || item1_p->P.GroupID == md->ActiveGroup))
                continue;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = item1_p->P.atomkind;

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
                        for (item2_p = md->SubDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(md->ActiveGroup == -1 || item2_p->P.GroupID == md->ActiveGroup))
                                continue;

                            if (item1_p != item2_p)
                            {
                                atomkind2 = item2_p->P.atomkind;
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

static void EAM_convert_r_to_r2(eam_t *eam, double *source, double *dest)
/* Consider two functions f1 and f2 with the relation f2(r^2)=f1(r)
 * between them. Given the array source[0..eam.Nr-1] containing a table
 * of the function f1, i.e. source[i]=f1(r_i), with r_i=i*eam.dr , this
 * routine returns an array dest[0..eam.Nr2-1] that contains the
 * values of the function f2 at points j*eam.dr2 . */
{
    double *sourceDD;
    int i;

    sourceDD = (double *)malloc(eam->Nr * sizeof(double));
    spline_prepare(eam->dr, source, eam->Nr, sourceDD);

    for (i=0; i < eam->Nr2; i++)
        dest[i] = spline_val(eam->dr, source, sourceDD, sqrt(i * eam->dr2));
    free(sourceDD);
}

static eam_t *load_DYNAMOsetfl(fmd_t *md, char *filePath)
{
    int i, j, k;

    eam_t *eam = (eam_t *)malloc(sizeof(eam_t));

    if (md->Is_MD_comm_root)
    {
        FILE *fp = fopen(filePath, "r");
        handleFileOpenError(fp, filePath);

        char str[1024];
        for (i=0; i<3; i++)
            fgets(str, 1024, fp);

        fscanf(fp, "%d", &eam->elementsNo);
        eam->elements = (eam_element_t *)malloc(eam->elementsNo * sizeof(eam_element_t));
        for (i=0; i < eam->elementsNo; i++)
            if (fscanf(fp, "%s", str) == 1)
            {
                eam->elements[i].name = (char *)malloc(strlen(str) + 1);
                strcpy(eam->elements[i].name, str);
            }

        double cutoff;
        fscanf(fp, "%d%lf%d%lf%lf", &eam->Nrho, &eam->drho, &eam->Nr, &eam->dr, &cutoff);
        eam->Nr2 = (eam->Nr += 2);
        assert( (eam->Nr-1) * eam->dr > cutoff );
        eam->cutoff_sqr = SQR(cutoff);
        eam->dr2 = SQR((eam->Nr-1) * eam->dr) / (eam->Nr2-1);

        double *tempArray = (double *)malloc(eam->Nr * sizeof(double));
        for (i=0; i < eam->elementsNo; i++)
        {
            eam->elements[i].eam = eam;
            fscanf(fp, "%s%lf%lf%s", str, &eam->elements[i].mass,
                &eam->elements[i].latticeParameter, str);
            eam->elements[i].mass /= MD_MASS_UNIT;
            eam->elements[i].F = (double *)malloc(eam->Nrho * sizeof(double));
            for (j=0; j < eam->Nrho; j++)
                fscanf(fp, "%lf", &eam->elements[i].F[j]);

            eam->elements[i].rho = (double *)malloc(eam->Nr2 * sizeof(double));
            for (j=0; j < eam->Nr-2; j++)  // read rho(r) values from file
                fscanf(fp, "%lf", &tempArray[j]);
            tempArray[eam->Nr-1] = tempArray[eam->Nr-2] = 0.;
            EAM_convert_r_to_r2(eam, tempArray, eam->elements[i].rho);

            eam->elements[i].phi = (double **)malloc(eam->elementsNo * sizeof(double *));
#ifdef USE_CSPLINE
            eam->elements[i].F_DD = (double *)malloc(eam->Nrho * sizeof(double));
            eam->elements[i].rhoDD = (double *)malloc(eam->Nr2 * sizeof(double));
            spline_prepare(eam->drho, eam->elements[i].F, eam->Nrho, eam->elements[i].F_DD);
            spline_prepare(eam->dr2, eam->elements[i].rho, eam->Nr2, eam->elements[i].rhoDD);
            eam->elements[i].phiDD = (double **)malloc(eam->elementsNo * sizeof(double *));
#endif
        }

        for (i=0; i < eam->elementsNo; i++)
            for (j=0; j<=i; j++)
            {
                for (k=0; k < eam->Nr-2; k++) // read r*phi values from file
                {
                    fscanf(fp, "%lf", &tempArray[k]);
                    if (k==0)
                        tempArray[k] = FLT_MAX;
                    else
                        tempArray[k] /= k * eam->dr;
                }
                tempArray[eam->Nr-1] = tempArray[eam->Nr-2] = 0.;
                eam->elements[i].phi[j] = (double *)malloc(eam->Nr2 * sizeof(double));
                EAM_convert_r_to_r2(eam, tempArray, eam->elements[i].phi[j]);
                eam->elements[j].phi[i] = eam->elements[i].phi[j];
#ifdef USE_CSPLINE
                eam->elements[i].phiDD[j] = (double *)malloc(eam->Nr2 * sizeof(double));
                spline_prepare(eam->dr2, eam->elements[i].phi[j], eam->Nr2, eam->elements[i].phiDD[j]);
                eam->elements[j].phiDD[i] = eam->elements[i].phiDD[j];
#endif
            }
        free(tempArray);
        fclose(fp);
    }

    MPI_Bcast(eam, sizeof(eam_t), MPI_CHAR, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
    if (md->SubDomain.myrank != ROOTPROCESS(md->SubDomain.numprocs))
    {
        eam->elements = (eam_element_t *)malloc(eam->elementsNo * sizeof(eam_element_t));
        for (i=0; i < eam->elementsNo; i++)
        {
            eam->elements[i].eam = eam;
            eam->elements[i].F = (double *)malloc(eam->Nrho * sizeof(double));
            eam->elements[i].rho = (double *)malloc(eam->Nr2 * sizeof(double));
            eam->elements[i].phi = (double **)malloc(eam->elementsNo * sizeof(double *));
            for (j=0; j<=i; j++)
            {
                eam->elements[i].phi[j] = (double *)malloc(eam->Nr2 * sizeof(double));
                eam->elements[j].phi[i] = eam->elements[i].phi[j];
            }
#ifdef USE_CSPLINE
            eam->elements[i].F_DD = (double *)malloc(eam->Nrho * sizeof(double));
            eam->elements[i].rhoDD = (double *)malloc(eam->Nr2 * sizeof(double));
            eam->elements[i].phiDD = (double **)malloc(eam->elementsNo * sizeof(double *));
            for (j=0; j<=i; j++)
            {
                eam->elements[i].phiDD[j] = (double *)malloc(eam->Nr2 * sizeof(double));
                eam->elements[j].phiDD[i] = eam->elements[i].phiDD[j];
            }
#endif
        }
    }

    for (i=0; i < eam->elementsNo; i++)
    {
        MPI_Bcast(&eam->elements[i].mass, 1, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        MPI_Bcast(&eam->elements[i].latticeParameter, 1, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

        unsigned namelen;
        if (md->Is_MD_comm_root)
            namelen = strlen(eam->elements[i].name);
        MPI_Bcast(&namelen, 1, MPI_UNSIGNED, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        if (md->SubDomain.myrank != ROOTPROCESS(md->SubDomain.numprocs))
            eam->elements[i].name = (char *)malloc(namelen + 1);
        MPI_Bcast(eam->elements[i].name, namelen+1, MPI_CHAR, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);

        MPI_Bcast(eam->elements[i].F, eam->Nrho, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        MPI_Bcast(eam->elements[i].rho, eam->Nr2, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        for (j=0; j<=i; j++)
            MPI_Bcast(eam->elements[i].phi[j], eam->Nr2, MPI_DOUBLE,
                      ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
#ifdef USE_CSPLINE
        MPI_Bcast(eam->elements[i].F_DD, eam->Nrho, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        MPI_Bcast(eam->elements[i].rhoDD, eam->Nr2, MPI_DOUBLE, ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
        for (j=0; j<=i; j++)
            MPI_Bcast(eam->elements[i].phiDD[j], eam->Nr2, MPI_DOUBLE,
                      ROOTPROCESS(md->SubDomain.numprocs), md->MD_comm);
#endif
    }

    return eam;
}

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_t *md, fmd_string_t filePath)
{
    eam_t *eam = load_DYNAMOsetfl(md, filePath);

    fmd_pot_t *pot = (fmd_pot_t *)malloc(sizeof(fmd_pot_t));
    pot->cat = POT_EAM_ALLOY;
    pot->data = eam;

    md->potsys.potlist = fmd_list_prepend(md->potsys.potlist, pot);

    return pot;
}

double fmd_pot_eam_getCutoffRadius(fmd_t *md, fmd_pot_t *pot)
{
    // TO-DO: handle error
    assert(pot->cat == POT_EAM_ALLOY);

    return sqrt(((eam_t *)pot->data)->cutoff_sqr);
}

void fmd_pot_eam_free(eam_t *eam)
{
    int i, j;

    for (i=0; i < eam->elementsNo; i++)
    {
        free(eam->elements[i].F);
        free(eam->elements[i].rho);
        for (j=0; j<=i; j++)
            free(eam->elements[i].phi[j]);
        free(eam->elements[i].phi);
#ifdef USE_CSPLINE
        free(eam->elements[i].F_DD);
        free(eam->elements[i].rhoDD);
        for (j=0; j<=i; j++)
            free(eam->elements[i].phiDD[j]);
        free(eam->elements[i].phiDD);
#endif
    }
    free(eam->elements);
    free(eam);
}

double fmd_pot_eam_getLatticeParameter(fmd_t *md, fmd_pot_t *pot, fmd_string_t element)
{
    // TO-DO: handle error
    assert(pot->cat == POT_EAM_ALLOY);

    eam_t *eam = (eam_t *)(pot->data);
    for (unsigned i=0; i < eam->elementsNo; i++)
        if (strcmp(element, eam->elements[i].name) == 0)
            return eam->elements[i].latticeParameter;

    // TO-DO: if element is not found in the potential, notify the library user
}

unsigned fmd_pot_eam_find_iloc(fmd_t *md, eam_t *eam, unsigned atomkind)
{
    for (int i=0; i < eam->elementsNo; i++)
        if (strcmp(md->potsys.atomkinds[atomkind].name, eam->elements[i].name) == 0)
            return i;

    // return -1 if the given eam potential doesn't include the specified atom kind
    return -1;
}
