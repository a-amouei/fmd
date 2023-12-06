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

#include <string.h>
#include <tgmath.h>
#include <float.h>
#include "eam.h"
#include "fmd-private.h"
#include "misc.h"
#include "list.h"
#include "general.h"

void _fmd_computeEAM_pass0(fmd_t *md, fmd_real_t FembSum)
{
    fmd_ituple_t jc, kc;
    int ir2, ir2_h;
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
    fmd_real_t PotEnergy = 0.0;

    /* iterate over all cells */
#pragma omp parallel for private(ic0,ic1,ic2,rho_i,rho_iDD,kc,jc,rv,r2,h,ir2, \
    ir2_h,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag) \
    shared(md,pottable) default(none) collapse(3) reduction(+:PotEnergy) schedule(static,1)
    for (ic0 = md->Subdomain.ic_start[0]; ic0 < md->Subdomain.ic_stop[0]; ic0++)
    for (ic1 = md->Subdomain.ic_start[1]; ic1 < md->Subdomain.ic_stop[1]; ic1++)
    for (ic2 = md->Subdomain.ic_start[2]; ic2 < md->Subdomain.ic_stop[2]; ic2++)
    {
        /* iterate over all particles in cell ic */

        cell_t *c1;
        unsigned i1;

        for (c1 = &md->Subdomain.grid[ic0][ic1][ic2], i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup)
                continue;

            for (int d=0; d<DIM; d++)
                FRC(c1, i1, d) = 0.0;

            eam_t *eam;
            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            /* iterate over neighbor cells of cell ic */
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0);
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1);
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2);

                        /* iterate over all particles in cell jc */

                        cell_t *c2;
                        unsigned i2;

                        for (c2 = &ARRAY_ELEMENT(md->Subdomain.grid, jc), i2=0; i2 < c2->parts_num; i2++)
                        {
                            if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup)
                                continue;

                            if ( (c1 != c2) || (i1 != i2) )
                            {
                                unsigned atomkind2 = c2->atomkind[i2];
                                fmd_real_t *x2 = &POS(c2, i2, 0);

                                EAM_PAIR_UPDATE_FORCE_AND_POTENERGY;
                            }
                        }
                    }
                }
            }
        }
    }

    PotEnergy = 0.5 * PotEnergy + FembSum;
    MPI_Allreduce(&PotEnergy, &md->GroupPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

void _fmd_computeEAM_pass1(fmd_t *md, fmd_real_t *FembSum_p)
{
    fmd_ituple_t jc, kc;
    int ir2, irho, ir2_h, irho_h;
    fmd_real_t r2;
    fmd_rtuple_t rv;
    fmd_real_t *rho, *rhoDD, *F, *F_DD;
    fmd_real_t a, b, h;
    int ic0, ic1, ic2;
    fmd_real_t Femb_sum=0;
    potpair_t **pottable = md->potsys.pottable;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    /* iterate over all cells */
    #pragma omp parallel for private(ic0,ic1,ic2,kc,jc,rv,r2,h,ir2,ir2_h,a,b,rho, \
      rhoDD,F,F_DD,irho,irho_h) shared(md,pottable,atomkinds) default(none) collapse(3) reduction(+:Femb_sum) \
      schedule(static,1)
    for (ic0 = md->Subdomain.ic_start[0]; ic0 < md->Subdomain.ic_stop[0]; ic0++)
    for (ic1 = md->Subdomain.ic_start[1]; ic1 < md->Subdomain.ic_stop[1]; ic1++)
    for (ic2 = md->Subdomain.ic_start[2]; ic2 < md->Subdomain.ic_stop[2]; ic2++)
    {
        /* iterate over all particles in cell ic */

        cell_t *c1;
        unsigned i1;

        for (c1 = &md->Subdomain.grid[ic0][ic1][ic2], i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup)
                continue;

            eam_t *eam;
            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);
            fmd_real_t rho_host = 0.0;

            /* iterate over neighbor cells of cell ic */
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0);
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1);
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2);

                        /* iterate over particles in cell jc */

                        cell_t *c2;
                        unsigned i2;

                        for (c2 = &ARRAY_ELEMENT(md->Subdomain.grid, jc), i2=0; i2 < c2->parts_num; i2++)
                        {
                            if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup)
                                continue;

                            if ( (c1 != c2) || (i1 != i2) )
                            {
                                unsigned atomkind2 = c2->atomkind[i2];
                                fmd_real_t *x2 = &POS(c2, i2, 0);

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

static void EAM_convert_r_to_r2(eam_t *eam, fmd_real_t *source, fmd_real_t *dest)
/* Consider two functions f1 and f2 with the relation f2(r^2)=f1(r)
 * between them. Given the array source[0..eam.Nr-1] containing a table
 * of the function f1, i.e. source[i]=f1(r_i), with r_i=i*eam.dr , this
 * routine returns an array dest[0..eam.Nr2-1] that contains the
 * values of the function f2 at points j*eam.dr2 . */
{
    fmd_real_t *SourceDD;
    int i;

    SourceDD = (fmd_real_t *)m_alloc(eam->Nr * sizeof(fmd_real_t));
    spline_prepare(eam->dr, source, eam->Nr, SourceDD);

    for (i=0; i < eam->Nr2; i++)
        dest[i] = spline_val(eam->dr, source, SourceDD, sqrt(i * eam->dr2));
    free(SourceDD);
}

static void create_mpi_eam(fmd_t *md, MPI_Datatype *mpi_eam)
{
    MPI_Datatype temptype;

    int eamblocklen[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    MPI_Aint eamdisplc[8] = {offsetof(eam_t, drho),
                             offsetof(eam_t, dr),
                             offsetof(eam_t, dr2),
                             offsetof(eam_t, cutoff_sqr),
                             offsetof(eam_t, ElementsNo),
                             offsetof(eam_t, Nrho),
                             offsetof(eam_t, Nr),
                             offsetof(eam_t, Nr2)};

    MPI_Datatype eamtype[8] = {FMD_MPI_REAL,
                               FMD_MPI_REAL,
                               FMD_MPI_REAL,
                               FMD_MPI_REAL,
                               MPI_INT,
                               MPI_INT,
                               MPI_INT,
                               MPI_INT};

    MPI_Type_create_struct(8, eamblocklen, eamdisplc, eamtype, &temptype);

    MPI_Type_create_resized(temptype, 0, sizeof(eam_t), mpi_eam);

    MPI_Type_free(&temptype);

    MPI_Type_commit(mpi_eam);
}

static eam_t *load_DYNAMOsetfl(fmd_t *md, char *FilePath)
{
    eam_t *eam = (eam_t *)m_alloc(sizeof(eam_t));

    if (md->Is_MD_comm_root)
    {
        FILE *fp = f_open(FilePath, "r");

        char str[1024];

        for (int i=0; i<3; i++)
            assert( fgets(str, 1024, fp) != NULL );

        assert( fscanf(fp, "%d", &eam->ElementsNo) == 1 );

        eam->elements = (eam_element_t *)m_alloc(eam->ElementsNo * sizeof(eam_element_t));

        for (int i=0; i < eam->ElementsNo; i++)
            if (fscanf(fp, "%s", str) == 1)
            {
                eam->elements[i].name = (char *)m_alloc(strlen(str) + 1);
                strcpy(eam->elements[i].name, str);
            }

        fmd_real_t cutoff;

        assert ( fscanf(fp, "%d%lf%d%lf%lf", &eam->Nrho, &eam->drho, &eam->Nr, &eam->dr, &cutoff) == 5 );
        eam->Nr2 = (eam->Nr += 2);
        assert( (eam->Nr-1) * eam->dr > cutoff );
        eam->cutoff_sqr = sqrr(cutoff);
        eam->dr2 = sqrr((eam->Nr-1) * eam->dr) / (eam->Nr2-1);

        fmd_real_t *TempArray = (fmd_real_t *)m_alloc(eam->Nr * sizeof(fmd_real_t));

        for (int i=0; i < eam->ElementsNo; i++)
        {
            eam->elements[i].eam = eam;
            assert ( fscanf(fp, "%s%lf%lf%s", str, &eam->elements[i].mass,
                &eam->elements[i].latticeParameter, str) == 4 );
            eam->elements[i].mass /= MD_MASS_UNIT;

            eam->elements[i].F = (fmd_real_t *)m_alloc(eam->Nrho * sizeof(fmd_real_t));
            for (int j=0; j < eam->Nrho; j++)
                assert( fscanf(fp, "%lf", &eam->elements[i].F[j]) == 1 );

            eam->elements[i].rho = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
            for (int j=0; j < eam->Nr-2; j++)  /* read rho(r) values from file */
                assert( fscanf(fp, "%lf", &TempArray[j]) == 1 );
            TempArray[eam->Nr-1] = TempArray[eam->Nr-2] = 0.;
            EAM_convert_r_to_r2(eam, TempArray, eam->elements[i].rho);

            eam->elements[i].phi = (fmd_real_t **)m_alloc(eam->ElementsNo * sizeof(fmd_real_t *));
#ifdef USE_CSPLINE
            eam->elements[i].F_DD = (fmd_real_t *)m_alloc(eam->Nrho * sizeof(fmd_real_t));
            eam->elements[i].rhoDD = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
            spline_prepare(eam->drho, eam->elements[i].F, eam->Nrho, eam->elements[i].F_DD);
            spline_prepare(eam->dr2, eam->elements[i].rho, eam->Nr2, eam->elements[i].rhoDD);
            eam->elements[i].phiDD = (fmd_real_t **)m_alloc(eam->ElementsNo * sizeof(fmd_real_t *));
#endif
        }

        for (int i=0; i < eam->ElementsNo; i++)
            for (int j=0; j<=i; j++)
            {
                for (int k=0; k < eam->Nr-2; k++) /* read r*phi values from file */
                {
                    assert( fscanf(fp, "%lf", &TempArray[k]) == 1 );

                    if (k==0)
                        TempArray[k] = FLT_MAX;
                    else
                        TempArray[k] /= k * eam->dr;
                }
                TempArray[eam->Nr-1] = TempArray[eam->Nr-2] = 0.;
                eam->elements[i].phi[j] = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
                EAM_convert_r_to_r2(eam, TempArray, eam->elements[i].phi[j]);
                eam->elements[j].phi[i] = eam->elements[i].phi[j];
#ifdef USE_CSPLINE
                eam->elements[i].phiDD[j] = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
                spline_prepare(eam->dr2, eam->elements[i].phi[j], eam->Nr2, eam->elements[i].phiDD[j]);
                eam->elements[j].phiDD[i] = eam->elements[i].phiDD[j];
#endif
            }

        free(TempArray);
        fclose(fp);
    }

    MPI_Datatype mpi_eam;

    create_mpi_eam(md, &mpi_eam);

    MPI_Bcast(eam, 1, mpi_eam, RANK0, md->MD_comm);

    MPI_Type_free(&mpi_eam);

    if (!md->Is_MD_comm_root)
    {
        eam->elements = (eam_element_t *)m_alloc(eam->ElementsNo * sizeof(eam_element_t));
        for (int i=0; i < eam->ElementsNo; i++)
        {
            eam->elements[i].eam = eam;
            eam->elements[i].F = (fmd_real_t *)m_alloc(eam->Nrho * sizeof(fmd_real_t));
            eam->elements[i].rho = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
            eam->elements[i].phi = (fmd_real_t **)m_alloc(eam->ElementsNo * sizeof(fmd_real_t *));

            for (int j=0; j<=i; j++)
            {
                eam->elements[i].phi[j] = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
                eam->elements[j].phi[i] = eam->elements[i].phi[j];
            }

#ifdef USE_CSPLINE
            eam->elements[i].F_DD = (fmd_real_t *)m_alloc(eam->Nrho * sizeof(fmd_real_t));
            eam->elements[i].rhoDD = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
            eam->elements[i].phiDD = (fmd_real_t **)m_alloc(eam->ElementsNo * sizeof(fmd_real_t *));
            for (int j=0; j<=i; j++)
            {
                eam->elements[i].phiDD[j] = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));
                eam->elements[j].phiDD[i] = eam->elements[i].phiDD[j];
            }
#endif
        }
    }

    for (int i=0; i < eam->ElementsNo; i++)
    {
        MPI_Bcast(&eam->elements[i].mass, 1, FMD_MPI_REAL, RANK0, md->MD_comm);
        MPI_Bcast(&eam->elements[i].latticeParameter, 1, FMD_MPI_REAL, RANK0, md->MD_comm);

        unsigned namelen;

        if (md->Is_MD_comm_root)
            namelen = strlen(eam->elements[i].name);

        MPI_Bcast(&namelen, 1, MPI_UNSIGNED, RANK0, md->MD_comm);

        if (!md->Is_MD_comm_root)
            eam->elements[i].name = (char *)m_alloc(namelen + 1);

        MPI_Bcast(eam->elements[i].name, namelen+1, MPI_CHAR, RANK0, md->MD_comm);

        MPI_Bcast(eam->elements[i].F, eam->Nrho, FMD_MPI_REAL, RANK0, md->MD_comm);
        MPI_Bcast(eam->elements[i].rho, eam->Nr2, FMD_MPI_REAL, RANK0, md->MD_comm);
        for (int j=0; j<=i; j++)
            MPI_Bcast(eam->elements[i].phi[j], eam->Nr2, FMD_MPI_REAL, RANK0, md->MD_comm);
#ifdef USE_CSPLINE
        MPI_Bcast(eam->elements[i].F_DD, eam->Nrho, FMD_MPI_REAL, RANK0, md->MD_comm);

        MPI_Bcast(eam->elements[i].rhoDD, eam->Nr2, FMD_MPI_REAL, RANK0, md->MD_comm);
        for (int j=0; j<=i; j++)
            MPI_Bcast(eam->elements[i].phiDD[j], eam->Nr2, FMD_MPI_REAL, RANK0, md->MD_comm);
#endif
    }

    return eam;
}

fmd_pot_t *fmd_pot_eam_alloy_load(fmd_t *md, fmd_string_t path)
{
    eam_t *eam = load_DYNAMOsetfl(md, path);

    fmd_pot_t *pot = (fmd_pot_t *)m_alloc(sizeof(fmd_pot_t));
    pot->cat = POT_EAM_ALLOY;
    pot->data = eam;

    md->potsys.potlist = fmd_list_prepend(md->potsys.potlist, pot);

    return pot;
}

fmd_real_t fmd_pot_eam_getCutoffRadius(fmd_t *md, fmd_pot_t *pot)
{
    // TO-DO: handle error
    assert(pot->cat == POT_EAM_ALLOY);

    return sqrt(((eam_t *)pot->data)->cutoff_sqr);
}

void _fmd_pot_eam_free(eam_t *eam)
{
    int i, j;

    for (i=0; i < eam->ElementsNo; i++)
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

fmd_real_t fmd_pot_eam_getLatticeParameter(fmd_t *md, fmd_pot_t *pot, fmd_string_t element)
{
    // TO-DO: handle error
    assert(pot->cat == POT_EAM_ALLOY);

    eam_t *eam = (eam_t *)(pot->data);
    for (unsigned i=0; i < eam->ElementsNo; i++)
        if (strcmp(element, eam->elements[i].name) == 0)
            return eam->elements[i].latticeParameter;

    // TO-DO: if element is not found in the potential, notify the library user
}

unsigned _fmd_pot_eam_find_iloc(fmd_t *md, eam_t *eam, unsigned atomkind)
{
    for (int i=0; i < eam->ElementsNo; i++)
        if (strcmp(md->potsys.atomkinds[atomkind].name, eam->elements[i].name) == 0)
            return i;

    // return -1 if the given eam potential doesn't include the specified atom kind
    return -1;
}
