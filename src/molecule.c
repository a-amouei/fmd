/*
  molecule.c: This file is part of Free Molecular Dynamics

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

#include "molecule.h"
#include "potential.h"
#include "base.h"
#include "array.h"
#include "list.h"

fmd_handle_t fmd_bond_addKind(fmd_t *md, fmd_bond_t cat, fmd_real_t coeffs[])
{
    int i = md->potsys.bondkinds_num;

    md->potsys.bondkinds = (bondkindp_t *)realloc(md->potsys.bondkinds,
      (i+1) * sizeof(bondkindp_t));

    switch (cat)
    {
        case FMD_BOND_HARMONIC:
            ;   // null statement
            bondkind_harmonic_t *bkh = (bondkind_harmonic_t *)malloc(sizeof(bondkind_harmonic_t));
            bkh->cat = cat;
            bkh->k = coeffs[0];
            bkh->r0 = coeffs[1];
            md->potsys.bondkinds[i] = (bondkindp_t)bkh;
            break;
    }

    md->potsys.bondkinds_num++;
    return i;
}

void fmd_bond_freeKinds(fmd_t *md)
{
    for (unsigned i=0; i < md->potsys.bondkinds_num; i++)
        free(md->potsys.bondkinds[i]);
    free(md->potsys.bondkinds);
    md->potsys.bondkinds_num = 0;
}

static int compare_LocalID_in_mkln(const void *a, const void *b)
{
    if ( ((molkind_atom_neighbor_t *)a)->atom->LocalID == *( (unsigned *)b ) )
        return 0;
    else
        return 1;
}

void fmd_bond_apply(fmd_t *md, fmd_handle_t bondkind, fmd_handle_t molkind,
  unsigned atom1, unsigned atom2)
{
    // TO-DO
    assert(atom1 != atom2);

    molkind_t *mk = &md->potsys.molkinds[molkind];
    list_t *ln1 = mk->atoms[atom1].neighbors;
    list_t *ln2 = mk->atoms[atom2].neighbors;

    // check one of the two 'neighbors' lists to see if it already contains the other atom
    list_t *item1 = fmd_list_find_custom(ln1, &atom2, compare_LocalID_in_mkln);
    list_t *item2;  // here in 'item1' & 'item2', '1' & '2' refer to list 1 & list 2

    if (item1 != NULL) // there is already a bond between the two atoms
    {
        // just find item2 as well and update the bond field of the
        // related molkind_atom_neighbor_t structures

        item2 = fmd_list_find_custom(ln2, &atom1, compare_LocalID_in_mkln);
        assert(item2 != NULL);

        ((molkind_atom_neighbor_t *)(item1->data))->bond = md->potsys.bondkinds[bondkind];
        ((molkind_atom_neighbor_t *)(item2->data))->bond = md->potsys.bondkinds[bondkind];
    }
    else // there is NOT already a bond between the two atoms
    {
        // create one!

        // here in 'data1' & 'data2', '1' & '2' refer to list 1 & list 2
        molkind_atom_neighbor_t *data1, *data2;

        // for neighbor list 1
        data1 = (molkind_atom_neighbor_t *)malloc(sizeof(molkind_atom_neighbor_t));
        data1->atom = &mk->atoms[atom2];
        data1->bond = md->potsys.bondkinds[bondkind];
        mk->atoms[atom1].neighbors = fmd_list_prepend(ln1, data1);

        // for neighbor list 2
        data2 = (molkind_atom_neighbor_t *)malloc(sizeof(molkind_atom_neighbor_t));
        data2->atom = &mk->atoms[atom1];
        data2->bond = md->potsys.bondkinds[bondkind];
        mk->atoms[atom2].neighbors = fmd_list_prepend(ln2, data2);
    }
}

void fmd_molecule_freeKinds()
{
    // TO-DO
    //free(md->potsys.molkinds + 1);  // the +1 shift is important here ( see fmd_molecule_addKind() )
}

fmd_handle_t fmd_molecule_addKind(fmd_t *md, fmd_string_t name, unsigned AtomsNum,
  unsigned AtomKinds[], fmd_rtuple_t AtomPositions[])
{
    int i = md->potsys.molkinds_num;

    // molkinds data start from index 1 because molkind=0 has special meaning
    molkind_t *genuine_pointer = (md->potsys.molkinds == NULL ? NULL : md->potsys.molkinds+1);
    md->potsys.molkinds = (molkind_t *)realloc(genuine_pointer, (i+1) * sizeof(molkind_t)) - 1;

    molkind_t *mk = &md->potsys.molkinds[i+1];
    mk->atoms_num = AtomsNum;
    size_t len = strlen(name);
    mk->name = (char *)malloc(len + 1);
    strcpy(mk->name, name);
    mk->distances = (unsigned **)_fmd_array_neat2d_create(AtomsNum, AtomsNum, sizeof(unsigned));

    mk->atoms = (molkind_atom_t *)malloc(AtomsNum * sizeof(molkind_atom_t));
    for (unsigned j=0; j<AtomsNum; j++)
    {
        mk->atoms[j].LocalID = j;
        mk->atoms[j].atomkind = AtomKinds[j];
        mk->atoms[j].neighbors = NULL;
        for (int d=0; d<3; d++)
            mk->atoms[j].position[d] = AtomPositions[j][d];
    }

    md->potsys.molkinds_num++;
    return i+1;
}

inline int find_neighbor_in_cell(cell_t *c, unsigned MolID, unsigned neighborID)
{
    for (int i=0; i < c->parts_num; i++)
        if (c->parts[i].core.MolID == MolID && c->parts[i].core.AtomID_local == neighborID)
            return i;

    return -1;
}

#define MAP_kc_TO_jc(kc, jc, d)                                             \
    if ((kc)[(d)] < md->SubDomain.ic_start[(d)])                            \
    {                                                                       \
        if (md->PBC[(d)] && md->ns[(d)] == 1)                               \
        {                                                                   \
            (jc)[(d)] = (kc)[(d)] + md->nc[(d)];                            \
            map_done=1;                                                     \
        }                                                                   \
        else                                                                \
            map_done=0;                                                     \
    }                                                                       \
    else if ((kc)[(d)] >= md->SubDomain.ic_stop[(d)])                       \
    {                                                                       \
        if (md->PBC[(d)] && md->ns[(d)] == 1)                               \
        {                                                                   \
            (jc)[(d)] = (kc)[(d)] - md->nc[(d)];                            \
            map_done=1;                                                     \
        }                                                                   \
        else                                                                \
            map_done=0;                                                     \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        (jc)[(d)] = (kc)[(d)];                                              \
        map_done=1;                                                         \
    }

#define MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, d)                                 \
    if ((kc)[(d)] < md->SubDomain.ic_start[(d)])                            \
    {                                                                       \
        if (md->PBC[(d)] && md->ns[(d)] == 1)                               \
            (jc)[(d)] = (kc)[(d)] + md->nc[(d)];                            \
        else                                                                \
            continue;                                                       \
    }                                                                       \
    else if ((kc)[(d)] >= md->SubDomain.ic_stop[(d)])                       \
    {                                                                       \
        if (md->PBC[(d)] && md->ns[(d)] == 1)                               \
            (jc)[(d)] = (kc)[(d)] - md->nc[(d)];                            \
        else                                                                \
            continue;                                                       \
    }                                                                       \
    else                                                                    \
        (jc)[(d)] = (kc)[(d)];

#define SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(xc)             \
    do                                                                      \
    {                                                                       \
        cell_t *cell = &md->SubDomain.grid[(xc)[0]][(xc)[1]][(xc)[2]];      \
        ind = find_neighbor_in_cell(cell, MolID, neighborID);               \
                                                                            \
        if (ind > -1) /* if the described neighbor was found */             \
        {                                                                   \
            *cnb = cell;                                                    \
            *index_nb = ind;                                                \
            return;                                                         \
        }                                                                   \
    } while (0)

/* the last two parameters are outputs of the function:
   (they determine the place of the neighbor particle in grid, if any found)

   (*c)[*index_nb] is the found neighbor.

   if index_nb equals -1, it means no neighbor is found                         */
static void find_neighbor(fmd_t *md, fmd_ituple_t ic, unsigned MolID, unsigned neighborID,
                          cell_t **cnb, int *index_nb)
{
    /* calculate maximum distance */
    int max_dist = 0;
    for (int d=0; d<3; d++)
    {
        int max_d;

        if (md->PBC[d] && md->ns[d] == 1)
            max_d = md->SubDomain.cell_num_nonmarg[d] / 2;
        else
        {
            max_d = ic[d] - md->SubDomain.ic_start[d];             // left
            unsigned tempo = md->SubDomain.ic_stop[d] - ic[d] - 1; // right
            if (tempo > max_d) max_d = tempo;
        }

        if (max_d > max_dist) max_dist = max_d;
    }

    int ind;

    /* treat dist=0 separately */
    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(ic);

    /* other dist values */
    for (int dist=1; dist <= max_dist; dist++)
    {
        fmd_ituple_t jc, kc;
        fmd_bool_t map_done;

        /* segment 1 out of 6 */
        kc[0] = ic[0] + dist;
        MAP_kc_TO_jc(kc, jc, 0);
        if (map_done)
            for (kc[1]=ic[1]-dist; kc[1]<=ic[1]+dist; kc[1]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(jc);
                }
            }

        /* segment 2 out of 6 */
        kc[0] = ic[0] - dist;
        MAP_kc_TO_jc(kc, jc, 0);
        if (map_done)
            for (kc[1]=ic[1]-dist; kc[1]<=ic[1]+dist; kc[1]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(jc);
                }
            }

        /* segment 3 out of 6 */
        kc[1] = ic[1] + dist;
        MAP_kc_TO_jc(kc, jc, 1);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(jc);
                }
            }

        /* segment 4 out of 6 */
        kc[1] = ic[1] - dist;
        MAP_kc_TO_jc(kc, jc, 1);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(jc);
                }
            }

        /* segment 5 out of 6 */
        kc[2] = ic[2] + dist;
        MAP_kc_TO_jc(kc, jc, 2);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[1]=ic[1]-(dist-1); kc[1]<=ic[1]+(dist-1); kc[1]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(jc);
                }
            }

        /* segment 6 out of 6 */
        kc[2] = ic[2] - dist;
        MAP_kc_TO_jc(kc, jc, 2);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[1]=ic[1]-(dist-1); kc[1]<=ic[1]+(dist-1); kc[1]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                    SEARCH_FOR_NEIGHBOR_IN_CELL_AND_SET_OUTPUT_IF_FOUND(jc);
                }
            }
    }

    *index_nb = -1;  /* if the specified neighbor atom was not found */
}

static int compare_LocalID_in_ln(const void *a, const void *b)
{
    if ( ((mol_atom_neighbor_t *)a)->LocalID == *( (unsigned *)b ) )
        return 0;
    else
        return 1;
}

void _fmd_matt_updateAtomNeighbors(fmd_t *md)
{
    fmd_ituple_t ic;
    cell_t *cell;
    int i;

    /* iterate over all particles in current subdomain */
    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (cell = &md->SubDomain.grid[ic[0]][ic[1]][ic[2]], i=0; i < cell->parts_num; i++)
        {
            particle_t *p = &cell->parts[i];

            if (p->core.molkind == 0) continue;

            /* iterate over all neighbor atoms in the template for the molecule */
            list_t *mkln = md->potsys.molkinds[p->core.molkind].atoms[p->core.AtomID_local].neighbors;

            while (mkln != NULL)
            {
                unsigned nblocal = ((molkind_atom_neighbor_t *)mkln->data)->atom->LocalID;

                /* does "neighbors" list of atom p contain an entry for a neighbor atom with local ID nblocal? */
                list_t *res = fmd_list_find_custom(p->neighbors, &nblocal, compare_LocalID_in_ln);

                if (res == NULL) /* NO */
                {
                    /* create that entry */
                    mol_atom_neighbor_t *man = (mol_atom_neighbor_t *)malloc(sizeof(mol_atom_neighbor_t));
                    man->LocalID = nblocal;
                    man->bond = ((molkind_atom_neighbor_t *)mkln->data)->bond;

                    /* the following two variables will describe the place of the neighbor atom */
                    cell_t *cnb;
                    int index;

                    find_neighbor(md, ic, p->core.MolID, nblocal, &cnb, &index);
                    man->cell = cnb;
                    man->index = index;
                    p->neighbors = fmd_list_prepend(p->neighbors, man);

                    /* does a neighbor atom with local ID nblocal actually exist in current subdomain? */
                    if (index != -1)  /* neighbor atom was found indeed */
                    {
                        /* also add an entry in the neighbor list of that neighbor atom and update it */
                        man = (mol_atom_neighbor_t *)malloc(sizeof(mol_atom_neighbor_t));
                        man->LocalID = p->core.AtomID_local;
                        man->bond = ((molkind_atom_neighbor_t *)mkln->data)->bond;
                        man->index = i;
                        man->cell = cell;
                        cnb->parts[index].neighbors = fmd_list_prepend(cnb->parts[index].neighbors, man);
                    }
                }
                else /* YES. It contains such entry. */
                {
                    /* does that entry refer to an atom existing in current subdomain? */

                    if ( ((mol_atom_neighbor_t *)res->data)->index == -1 ) /* no, it refers to nothing. */
                    {
                        /* try to find that neighbor atom in the subdomain, and update the entry accordingly */

                        /* the following two variables will describe the place of the neighbor atom */
                        cell_t *cnb;
                        int index;

                        find_neighbor(md, ic, p->core.MolID, nblocal, &cnb, &index);
                        ((mol_atom_neighbor_t *)res->data)->cell = cnb;
                        ((mol_atom_neighbor_t *)res->data)->index = index;

                        if (index != -1)  /* neighbor atom was found */
                        {
                            /* also update the corresponding entry in the neighbors list of the neighbor atom */
                            list_t *f = fmd_list_find_custom(cnb->parts[index].neighbors, &p->core.AtomID_local,
                                                             compare_LocalID_in_ln);
                            ((mol_atom_neighbor_t *)f->data)->cell = cell;
                            ((mol_atom_neighbor_t *)f->data)->index = i;
                        }
                    }
                }

                mkln = mkln->next;
            }
        }
}

#define _COMPUTE_rv_AND_r2(p1, p2, r0)                                       \
    for (d=0; d<3; d++)                                                      \
    {                                                                        \
        rv[d] = (p1)->core.x[d] - (p2)->core.x[d];                           \
        if (md->ns[d] == 1)                                                  \
        {                                                                    \
            if (rv[d] > 2*(r0))                                              \
                rv[d] -= md->l[d];                                           \
            else if (rv[d] < -2*(r0))                                        \
                rv[d] += md->l[d];                                           \
        }                                                                    \
    }                                                                        \
    r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);

void fmd_dync_computeBondForce(fmd_t *md)
{
    int ic0, ic1, ic2;
    fmd_real_t PotEnergy = 0.0;

    /* iterate over all cells(lists) */
    #pragma omp parallel for private(ic0,ic1,ic2) \
      shared(md) default(none) collapse(3) reduction(+:PotEnergy) schedule(static,1)
    for (ic0 = md->SubDomain.ic_start[0]; ic0 < md->SubDomain.ic_stop[0]; ic0++)
        for (ic1 = md->SubDomain.ic_start[1]; ic1 < md->SubDomain.ic_stop[1]; ic1++)
            for (ic2 = md->SubDomain.ic_start[2]; ic2 < md->SubDomain.ic_stop[2]; ic2++)
            {
                cell_t *cell;
                int i;

                /* iterate over all atoms in cell ic */
                for (cell = &md->SubDomain.grid[ic0][ic1][ic2], i=0; i < cell->parts_num; i++)
                {
                    particle_t *p1 = &cell->parts[i];

                    if (p1->core.molkind != 0)
                    {
                        int index = ((mol_atom_neighbor_t *)(p1->neighbors->data))->index;
                        particle_t *p2 = &((mol_atom_neighbor_t *)(p1->neighbors->data))->cell->parts[index];;

                        if (p1->core.AtomID_local > p2->core.AtomID_local)
                        {
                            fmd_real_t r2;
                            fmd_rtuple_t rv;
                            int d;

                            bondkind_harmonic_t *bond = (bondkind_harmonic_t *)((mol_atom_neighbor_t *)(p1->neighbors->data))->bond;
                            fmd_real_t r0 = bond->r0;
                            _COMPUTE_rv_AND_r2(p1, p2, r0);
                            fmd_real_t k = bond->k;
                            fmd_real_t r = sqrt(r2);
                            fmd_rtuple_t vek;
                            for (d=0; d<3; d++)
                            {
                                vek[d] = -2*k*(r-r0)*rv[d]/r;
                                p1->F[d] += vek[d];
                                p2->F[d] -= vek[d];
                            }
                            PotEnergy += k*(r-r0)*(r-r0);
                        }
                    }
                }
            }
    md->TotalPotentialEnergy += PotEnergy;
}
