/*
  potential.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

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

// all functions for loading potential files and preparing or managing data structures related to
// potentials come here.

#include "potential.h"
#include "base.h"
#include "array.h"
#include "list.h"
#include "eam.h"

void fmd_pot_setCutoffRadius(fmd_t *md, double cutoff)
{
    md->cutoffRadius = cutoff;
}

void fmd_matt_setAtomKinds(fmd_t *md, unsigned number, const fmd_string_t names[], const double masses[])
{
    // TO-DO: error should be handled here
    assert(number > 0);

    md->potsys.atomkinds_num = number;
    md->potsys.atomkinds = (atomkind_t *)malloc(number * sizeof(atomkind_t));
    for (unsigned i=0; i<number; i++)
    {
        md->potsys.atomkinds[i].mass = masses[i] / MD_MASS_UNIT;  // convert from amu to internal mass unit
        size_t len = strlen(names[i]);
        md->potsys.atomkinds[i].name = (char *)malloc(len + 1);
        strcpy(md->potsys.atomkinds[i].name, names[i]);
    }
}

static void potlist_free(fmd_t *md)
{
    list_t *potlist = md->potsys.potlist;

    while (potlist != NULL)
    {
        fmd_pot_t *pot = (fmd_pot_t *)(potlist->data);

        switch (pot->kind)
        {
            case POTKIND_EAM_ALLOY:
                fmd_pot_eam_free((eam_t *)(pot->data));
                break;

            default:
                free(pot->data);
        }

        potlist = potlist->next;
    }

    fmd_list_free(md->potsys.potlist);
}


void fmd_potsys_free(fmd_t *md)
{
    if (md->potsys.atomkinds != NULL)
    {
        free(md->potsys.atomkinds);
        md->potsys.atomkinds = NULL;
    }

    if (md->potsys.bonds != NULL)
    {
        free(md->potsys.bonds);
        md->potsys.bonds = NULL;
    }

    if (md->potsys.pottable != NULL)
    {
        fmd_array_neat2d_free((void **)md->potsys.pottable);
        md->potsys.pottable = NULL;
    }

    if (md->potsys.potlist != NULL)
    {
        potlist_free(md);
        md->potsys.potlist = NULL;
    }

    if (md->potsys.potkinds != NULL)
    {
        fmd_list_free(md->potsys.potkinds);
        md->potsys.potkinds = NULL;
    }
}

void fmd_potsys_init(fmd_t *md)
{
    md->potsys.atomkinds = NULL;
    md->potsys.potlist = NULL;
    md->potsys.pottable = NULL;
    md->potsys.potkinds = NULL;
    md->potsys.bonds = NULL;
}

static void pottable_create(fmd_t *md)
{
    md->potsys.pottable = (potpair_t **)fmd_array_neat2d_create(md->potsys.atomkinds_num,
                                                                  md->potsys.atomkinds_num,
                                                                  sizeof(potpair_t));
    for (unsigned i=0; i < md->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
            md->potsys.pottable[i][j].kind = POTKIND_NONE;
}

void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot)
{
    // create the pottable if doesn't exist
    if (md->potsys.pottable == NULL) pottable_create(md);

    if (pot->kind == POTKIND_EAM_ALLOY)
    {
        // find the local indices
        unsigned loc1, loc2;

        loc1 = fmd_pot_eam_find_iloc(md, pot->data, atomkind1);
        loc2 = fmd_pot_eam_find_iloc(md, pot->data, atomkind2);

        if (loc1 == -1 || loc2 == -1)
        {
            // TO-DO: error should be handled here
            assert(loc1 != -1 && loc2 != -1);
        }

        md->potsys.pottable[atomkind1][atomkind2].iloc = loc1;
        md->potsys.pottable[atomkind1][atomkind2].jloc = loc2;
        md->potsys.pottable[atomkind2][atomkind1].iloc = loc2;
        md->potsys.pottable[atomkind2][atomkind1].jloc = loc1;

        // set eam_element
        md->potsys.atomkinds[atomkind1].eam_element = &((eam_t *)pot->data)->elements[loc1];
        md->potsys.atomkinds[atomkind2].eam_element = &((eam_t *)pot->data)->elements[loc2];
    }

    //
    md->potsys.pottable[atomkind1][atomkind2].kind =
      md->potsys.pottable[atomkind2][atomkind1].kind = pot->kind;
    md->potsys.pottable[atomkind1][atomkind2].data =
      md->potsys.pottable[atomkind2][atomkind1].data = pot->data;
}

static int potkind_compare(const void *a, const void *b)
{
    if ( *( (potkind_t *)a ) == *( (potkind_t *)b ) )
        return 0;
    else
        return 1;
}

static void pot_hybridpasses_update(fmd_t *md)
{
    // fill hybridpasses with zeros
    memset(md->potsys.hybridpasses, 0, sizeof(md->potsys.hybridpasses));

    list_t *potkinds = md->potsys.potkinds;
    while (potkinds != NULL)
    {
        switch (*(potkind_t *)(potkinds->data))
        {
            case POTKIND_MORSE:
                md->potsys.hybridpasses[0] = 1;
                break;

            case POTKIND_LJ_6_12:
                md->potsys.hybridpasses[0] = 1;
                break;

            case POTKIND_EAM_ALLOY:
                md->potsys.hybridpasses[0] = 1;
                md->potsys.hybridpasses[1] = 1;
                break;
        }

        potkinds = potkinds->next;
    }
}

// TO-DO?: first, clean the potkinds list
void fmd_pot_prepareForForceComp(fmd_t *md)
{
    // TO-DO: error should be handled here
    assert(md->potsys.pottable != NULL);

    md->potsys.potkinds_num = 0;

    for (unsigned i=0; i < md->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
        {
            potpair_t *potpair = &md->potsys.pottable[i][j];

            // TO-DO: error should be handled here
            assert(potpair->kind != POTKIND_NONE);

            // add the potkind to potkinds list, if isn't already included there
            if (fmd_list_find_custom(md->potsys.potkinds, &potpair->kind, potkind_compare) == NULL)
            {
                md->potsys.potkinds_num++;
                md->potsys.potkinds = fmd_list_prepend(md->potsys.potkinds, &potpair->kind);
            }

            // do these two atomkinds use EAM?
            if (potpair->kind != POTKIND_EAM_ALLOY)
                md->potsys.atomkinds[i].eam_element =
                md->potsys.atomkinds[j].eam_element = NULL;
        }

    pot_hybridpasses_update(md);
}
