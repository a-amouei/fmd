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
    potpair_t **pottable = md->potsys.pottable;
    fmd_real_t PotEnergy = 0.0;

    _fmd_clean_forces(md);

    /* iterate over all non-margin cells */

    #pragma omp parallel for shared(md,pottable) default(none) reduction(+:PotEnergy) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic=0; ic < md->subd.nc; ic++)
    {
        cell_t *c1 = md->subd.grid + ic;

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            /* iterate over non-margin neighbor cells of cell c1 */

            for (int jc=0; jc < c1->cnb2len; jc++)
            {
                cell_t *c2 = c1->cnb2[jc];

                /* iterate over all particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    unsigned atomkind2 = c2->atomkind[i2];

                    EAM_PAIR_UPDATE_FORCE_AND_POTENERGY2(x1, atomkind1, c1, i1, atomkind2,
                                                         c2, i2, PotEnergy, pottable);
                }
            }

            cell_t *c2 = c1;

            /* iterate over particles in cell c2=c1 with i2 < i1 */

            for (int i2=0; i2 < i1; i2++)
            {
                if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                unsigned atomkind2 = c2->atomkind[i2];

                EAM_PAIR_UPDATE_FORCE_AND_POTENERGY2(x1, atomkind1, c1, i1, atomkind2,
                                                     c2, i2, PotEnergy, pottable);
            }
        }
    }

    /* iterate over all margin cells */

    #pragma omp parallel for shared(md,pottable) default(none) reduction(+:PotEnergy) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic=md->subd.nc; ic < md->subd.ncm; ic++)
    {
        cell_t *c1 = md->subd.grid + ic;

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            /* iterate over non-margin neighbor cells of cell c1 */

            for (int jc=0; jc < c1->cnb1len; jc++)
            {
                cell_t *c2 = c1->cnb1[jc];

                /* iterate over particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    unsigned atomkind2 = c2->atomkind[i2];

                    EAM_PAIR_UPDATE_FORCE(x1, atomkind1, c1, i1, atomkind2, c2, i2, pottable);
                }
            }

            for (int jc=0; jc < c1->cnb2len; jc++)
            {
                cell_t *c2 = c1->cnb2[jc];

                /* iterate over particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    unsigned atomkind2 = c2->atomkind[i2];

                    EAM_PAIR_UPDATE_FORCE_AND_POTENERGY1(x1, atomkind1, c1, i1, atomkind2,
                                                         c2, i2, PotEnergy, pottable);
                }
            }
        }
    }

    PotEnergy += FembSum;
    MPI_Allreduce(&PotEnergy, &md->GroupPotentialEnergy, 1, FMD_MPI_REAL, MPI_SUM, md->MD_comm);
}

void _fmd_clean_vaream(fmd_t *md)
{
    int i;
    cell_t *c;

    for (int ic=0; ic < md->subd.nc; ic++)
        for (i=0, c = md->subd.grid + ic; i < c->parts_num; i++)
            c->vaream[i] = 0.;

}

#ifdef USE_CSPLINE
#define EAM_COMPUTE_FembPrime_AND_UPDATE_FembSum(c1, i1, atomkind1, atomkinds, FembSum)    \
    do                                                                                     \
    {                                                                                      \
        eam_element_t *el = atomkinds[atomkind1].eam_element;                              \
        eam_t *eam = el->eam;                                                              \
        fmd_real_t h = eam->drho;                                                          \
        fmd_real_t rho_host = c1->vaream[i1];                                              \
        int irho = (int)(rho_host / h);                                                    \
        assert(irho < eam->Nrho - 1);                                                      \
        int irho_h = irho + 1;                                                             \
        fmd_real_t *F = el->F;                                                             \
        fmd_real_t *F_DD = el->F_DD;                                                       \
        fmd_real_t a = irho_h - rho_host/h;                                                \
        fmd_real_t b = 1-a;                                                                \
        c1->vaream[i1] = SPLINE_DERIV(a,b,F,irho,irho_h,F_DD,h);                           \
        FembSum += SPLINE_VAL(a,b,F,irho,irho_h,F_DD,h);                                   \
    } while (0)
#else
#define EAM_COMPUTE_FembPrime_AND_UPDATE_FembSum(c1, i1, atomkind1, atomkinds, FembSum)    \
    do                                                                                     \
    {                                                                                      \
        eam_element_t *el = atomkinds[atomkind1].eam_element;                              \
        eam_t *eam = el->eam;                                                              \
        fmd_real_t h = eam->drho;                                                          \
        fmd_real_t rho_host = c1->vaream[i1];                                              \
        int irho = (int)(rho_host / h);                                                    \
        assert(irho < eam->Nrho - 1);                                                      \
        int irho_h = irho + 1;                                                             \
        fmd_real_t *F = el->F;                                                             \
        c1->vaream[i1] = (F[irho_h] - F[irho]) / h;                                        \
        FembSum += F[irho] + (rho_host - irho * h) * c1->vaream[i1];                       \
    } while (0)
#endif

fmd_real_t _fmd_calcFembPrime(fmd_t *md)
{
    fmd_real_t FembSum = 0.;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    for (int ic=0; ic < md->subd.nc; ic++)
    {
        cell_t *c1 = md->subd.grid + ic;

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];

            if (atomkinds[atomkind1].eam_element != NULL)
                EAM_COMPUTE_FembPrime_AND_UPDATE_FembSum(c1, i1, atomkind1, atomkinds, FembSum);
        }
    }

    return FembSum;
}

void _fmd_computeEAM_pass1(fmd_t *md, fmd_real_t *FembSum)
{
    potpair_t **pottable = md->potsys.pottable;
    atomkind_t *atomkinds = md->potsys.atomkinds;

    _fmd_clean_vaream(md);

    /* iterate over all non-margin cells */

    #pragma omp parallel for shared(md,pottable,atomkinds) default(none) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic=0; ic < md->subd.nc; ic++)
    {
        cell_t *c1 = md->subd.grid + ic;

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            /* iterate over non-margin neighbor cells of cell c1 */

            for (int jc=0; jc < c1->cnb2len; jc++)
            {
                cell_t *c2 = c1->cnb2[jc];

                /* iterate over particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    unsigned atomkind2 = c2->atomkind[i2];

                    EAM_PAIR_UPDATE_rho_host2(x1, atomkind1, c1, i1, atomkind2, c2, i2, pottable);
                }
            }

            cell_t *c2 = c1;

            /* iterate over particles in cell c2=c1 with i2 < i1 */

            for (int i2=0; i2 < i1; i2++)
            {
                if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                unsigned atomkind2 = c2->atomkind[i2];

                EAM_PAIR_UPDATE_rho_host2(x1, atomkind1, c1, i1, atomkind2, c2, i2, pottable);
            }
        }
    }

    /* iterate over all margin cells */

    #pragma omp parallel for shared(md,pottable,atomkinds) default(none) \
      schedule(dynamic,1) num_threads(md->numthreads)

    for (int ic=md->subd.nc; ic < md->subd.ncm; ic++)
    {
        cell_t *c1 = md->subd.grid + ic;

        /* iterate over all particles in cell c1 */

        for (int i1=0; i1 < c1->parts_num; i1++)
        {
            if (md->ActiveGroup != FMD_GROUP_ALL && c1->GroupID[i1] != md->ActiveGroup) continue;

            unsigned atomkind1 = c1->atomkind[i1];
            fmd_real_t *x1 = &POS(c1, i1, 0);

            /* iterate over non-margin neighbor cells of cell c1 */

            for (int jc=0; jc < c1->cnb0len; jc++)
            {
                cell_t *c2 = c1->cnb0[jc];

                /* iterate over particles in cell c2 */

                for (int i2=0; i2 < c2->parts_num; i2++)
                {
                    if (md->ActiveGroup != FMD_GROUP_ALL && c2->GroupID[i2] != md->ActiveGroup) continue;

                    unsigned atomkind2 = c2->atomkind[i2];

                    EAM_PAIR_UPDATE_rho_host1(x1, atomkind1, c1, i1, atomkind2, c2, i2, pottable);
                }
            }
        }
    }

    *FembSum = _fmd_calcFembPrime(md);
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
        {
            char *ret = fgets(str, 1024, fp);
            assert(ret != NULL); /* TO-DO: handle error */
        }

        int numread = fscanf(fp, "%d", &eam->ElementsNo);
        assert(numread == 1); /* TO-DO: handle error */

        eam->elements = (eam_element_t *)m_alloc(eam->ElementsNo * sizeof(eam_element_t));

        for (int i=0; i < eam->ElementsNo; i++)
        {
            numread = fscanf(fp, "%s", str);
            assert(numread == 1); /* TO-DO: handle error */

            eam->elements[i].name = (char *)m_alloc(strlen(str) + 1);
            strcpy(eam->elements[i].name, str);
        }

        fmd_real_t cutoff;

        numread = fscanf(fp, "%d%lf%d%lf%lf", &eam->Nrho, &eam->drho, &eam->Nr, &eam->dr, &cutoff);
        assert(numread == 5); /* TO-DO: handle error */

        eam->Nr2 = (eam->Nr += 2);
        assert( (eam->Nr-1) * eam->dr > cutoff ); /* TO-DO: handle error */
        eam->cutoff_sqr = sqrr(cutoff);
        eam->dr2 = sqrr((eam->Nr-1) * eam->dr) / (eam->Nr2-1);

        fmd_real_t *TempArray = (fmd_real_t *)m_alloc(eam->Nr * sizeof(fmd_real_t));

        for (int i=0; i < eam->ElementsNo; i++)
        {
            eam->elements[i].eam = eam;

            numread = fscanf(fp, "%s%lf%lf%s", str, &eam->elements[i].mass,
                             &eam->elements[i].latticeParameter, str);
            assert(numread == 4); /* TO-DO: handle error */

            eam->elements[i].mass /= MD_MASS_UNIT;

            eam->elements[i].F = (fmd_real_t *)m_alloc(eam->Nrho * sizeof(fmd_real_t));

            for (int j=0; j < eam->Nrho; j++)
            {
                numread = fscanf(fp, "%lf", &eam->elements[i].F[j]);
                assert(numread == 1); /* TO-DO: handle error */
            }

            eam->elements[i].rho = (fmd_real_t *)m_alloc(eam->Nr2 * sizeof(fmd_real_t));

            for (int j=0; j < eam->Nr-2; j++)  /* read rho(r) values from file */
            {
                numread = fscanf(fp, "%lf", &TempArray[j]);
                assert(numread == 1); /* TO-DO: handle error */
            }
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
                    numread = fscanf(fp, "%lf", &TempArray[k]);
                    assert(numread == 1); /* TO-DO: handle error */

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

    md->potsys.potlist = _fmd_list_prepend(md->potsys.potlist, pot);

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
