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

unsigned fmd_bond_addKind(fmd_t *md, fmd_bond_t cat, double coeffs[])
{
    unsigned i = md->potsys.bondkinds_num;

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

static int LocalID_compare(const void *a, const void *b)
{
    if ( ((molkind_atom_neighbor_t *)a)->atom->LocalID == *( (unsigned *)b ) )
        return 0;
    else
        return 1;
}

void fmd_bond_apply(fmd_t *md, unsigned bondkind, unsigned molkind,
  unsigned atom1, unsigned atom2)
{
    // TO-DO
    assert(atom1 != atom2);

    molkind_t *mk = &md->potsys.molkinds[molkind];
    list_t *ln1 = mk->atoms[atom1].neighbors;
    list_t *ln2 = mk->atoms[atom2].neighbors;

    // check one of the two 'neighbors' lists to see if it already contains the other atom
    list_t *item1 = fmd_list_find_custom(ln1, &atom2, LocalID_compare);
    list_t *item2;  // here in 'item1' & 'item2', '1' & '2' refer to list 1 & list 2

    if (item1 != NULL) // there is already a bond between the two atoms
    {
        // just find item2 as well and update the bond field of the
        // related molkind_atom_neighbor_t structures

        item2 = fmd_list_find_custom(ln2, &atom1, LocalID_compare);
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

unsigned fmd_molecule_addKind(fmd_t *md, fmd_string_t name, unsigned AtomsNum,
  unsigned AtomKinds[], double AtomPositions[][3])
{
    unsigned i = md->potsys.molkinds_num;

    // molkinds data start from index 1 because molkind=0 has special meaning
    molkind_t *genuine_pointer = (md->potsys.molkinds == NULL ? NULL : md->potsys.molkinds+1);
    md->potsys.molkinds = (molkind_t *)realloc(genuine_pointer, (i+1) * sizeof(molkind_t)) - 1;

    molkind_t *mk = &md->potsys.molkinds[i+1];
    mk->atoms_num = AtomsNum;
    size_t len = strlen(name);
    mk->name = (char *)malloc(len + 1);
    strcpy(mk->name, name);
    mk->distances = (unsigned **)fmd_array_neat2d_create(AtomsNum, AtomsNum, sizeof(unsigned));

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

#define FIND_NEIGHBOR_IN_jc(jc)                                                                             \
    for (item_p = md->SubDomain.grid[(jc)[0]][(jc)[1]][(jc)[2]]; item_p != NULL; item_p = item_p->next_p)   \
        if (item_p->P.MolID == MolID && item_p->P.AtomID_local == neighborID)                               \
            break;

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

static TParticleListItem *find_neighbor(fmd_t *md, int ic[3], unsigned MolID, unsigned neighborID)
{
    // calculate maximum distance
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

    TParticleListItem *item_p;

    // treat dist=0 separately
    FIND_NEIGHBOR_IN_jc(ic);
    if (item_p != NULL) return item_p;

    // other dist values
    for (int dist=1; dist <= max_dist; dist++)
    {
        int jc[3], kc[3];
        fmd_bool_t map_done;

        // segment 1 out of 6
        kc[0] = ic[0] + dist;
        MAP_kc_TO_jc(kc, jc, 0);
        if (map_done)
            for (kc[1]=ic[1]-dist; kc[1]<=ic[1]+dist; kc[1]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    FIND_NEIGHBOR_IN_jc(jc);
                    if (item_p != NULL) return item_p;
                }
            }

        // segment 2 out of 6
        kc[0] = ic[0] - dist;
        MAP_kc_TO_jc(kc, jc, 0);
        if (map_done)
            for (kc[1]=ic[1]-dist; kc[1]<=ic[1]+dist; kc[1]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    FIND_NEIGHBOR_IN_jc(jc);
                    if (item_p != NULL) return item_p;
                }
            }

        // segment 3 out of 6
        kc[1] = ic[1] + dist;
        MAP_kc_TO_jc(kc, jc, 1);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    FIND_NEIGHBOR_IN_jc(jc);
                    if (item_p != NULL) return item_p;
                }
            }

        // segment 4 out of 6
        kc[1] = ic[1] - dist;
        MAP_kc_TO_jc(kc, jc, 1);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[2]=ic[2]-dist; kc[2]<=ic[2]+dist; kc[2]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 2);
                    FIND_NEIGHBOR_IN_jc(jc);
                    if (item_p != NULL) return item_p;
                }
            }

        // segment 5 out of 6
        kc[2] = ic[2] + dist;
        MAP_kc_TO_jc(kc, jc, 2);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[1]=ic[1]-(dist-1); kc[1]<=ic[1]+(dist-1); kc[1]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                    FIND_NEIGHBOR_IN_jc(jc);
                    if (item_p != NULL) return item_p;
                }
            }

        // segment 6 out of 6
        kc[2] = ic[2] - dist;
        MAP_kc_TO_jc(kc, jc, 2);
        if (map_done)
            for (kc[0]=ic[0]-(dist-1); kc[0]<=ic[0]+(dist-1); kc[0]++)
            {
                MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 0);
                for (kc[1]=ic[1]-(dist-1); kc[1]<=ic[1]+(dist-1); kc[1]++)
                {
                    MAP_kc_TO_jc_INSIDE_LOOP(kc, jc, 1);
                    FIND_NEIGHBOR_IN_jc(jc);
                    if (item_p != NULL) return item_p;
                }
            }
    }

    return NULL;
}

void fmd_matt_updateNeighbors(fmd_t *md)
{
    int ic[3];
    TParticleListItem *item_p;

    ITERATE(ic, md->SubDomain.ic_start, md->SubDomain.ic_stop)
        for (item_p = md->SubDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            list_t *mkln = md->potsys.molkinds[item_p->P.molkind].atoms[item_p->P.AtomID_local].neighbors;

            while (mkln != NULL)
            {
                unsigned nblocal = ((molkind_atom_neighbor_t *)mkln->data)->atom->LocalID;

                find_neighbor(md, ic, item_p->P.MolID, nblocal);

                mkln = mkln->next;
            }
        }
}
