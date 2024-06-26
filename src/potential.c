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

#include <string.h>
#include <tgmath.h>
#include "potential.h"
#include "fmd-private.h"
#include "misc.h"
#include "array.h"
#include "list.h"
#include "eam.h"
#include "molecule.h"
#include "general.h"
#include "morse.h"
#include "lj.h"

void fmd_matt_setAtomKinds(fmd_t *md, unsigned number, const fmd_string_t names[], const fmd_real_t masses[])
{
    // TO-DO: error should be handled here
    assert(number > 0);

    md->potsys.atomkinds_num = number;
    md->potsys.atomkinds = (atomkind_t *)m_alloc(number * sizeof(atomkind_t));

    for (unsigned i=0; i<number; i++)
    {
        md->potsys.atomkinds[i].mass = masses[i] / MD_MASS_UNIT;  // convert from amu to internal mass unit
        size_t len = strlen(names[i]);
        md->potsys.atomkinds[i].name = (char *)m_alloc(len + 1);
        strcpy(md->potsys.atomkinds[i].name, names[i]);

        md->potsys.atomkinds[i].eam_element = NULL;
    }
}

static void potlist_free(fmd_t *md)
{
    list_t *potlist = md->potsys.potlist;

    while (potlist != NULL)
    {
        fmd_pot_t *pot = (fmd_pot_t *)(potlist->data);

        switch (pot->cat)
        {
            case POT_EAM_ALLOY:
                _fmd_pot_eam_free((eam_t *)(pot->data));
                break;

            default:
                free(pot->data);
        }

        potlist = potlist->next;
    }

    _fmd_list_free(md->potsys.potlist);
    md->potsys.potlist = NULL;
}

static void atomkinds_free(fmd_t *md)
{
    for (int i=0; i < md->potsys.atomkinds_num; i++)
        free(md->potsys.atomkinds[i].name);

    free(md->potsys.atomkinds);
    md->potsys.atomkinds_num = 0;
}

void _fmd_potsys_free(fmd_t *md)
{
    if (md->potsys.molkinds != NULL)
    {
        // TO-DO: call fmd_molecule_freeKinds()
        md->potsys.molkinds = NULL;
    }

    if (md->potsys.atomkinds != NULL)
    {
        atomkinds_free(md);
        md->potsys.atomkinds = NULL;
    }

    if (md->potsys.bondkinds != NULL)
    {
        fmd_bond_freeKinds(md);
        md->potsys.bondkinds = NULL;
    }

    if (md->potsys.pottable != NULL)
    {
        _fmd_array_neat2d_free((void **)md->potsys.pottable);
        md->potsys.pottable = NULL;
    }

    if (md->potsys.potlist != NULL)
    {
        potlist_free(md);
        md->potsys.potlist = NULL;
    }

    if (md->potsys.potcats != NULL)
    {
        _fmd_list_free(md->potsys.potcats);
        md->potsys.potcats = NULL;
    }
}

void _fmd_potsys_init(fmd_t *md)
{
    md->potsys.atomkinds = NULL;
    md->potsys.potlist = NULL;
    md->potsys.pottable = NULL;
    md->potsys.potcats = NULL;
    md->potsys.bondkinds = NULL;
    md->potsys.bondkinds_num = 0;
    md->potsys.molkinds = NULL;
    md->potsys.molkinds_num = 0;
}

static void pottable_create(fmd_t *md)
{
    md->potsys.pottable = (potpair_t **)_fmd_array_neat2d_create(md->potsys.atomkinds_num,
                                                                 md->potsys.atomkinds_num,
                                                                 sizeof(potpair_t));
    for (unsigned i=0; i < md->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
            md->potsys.pottable[i][j].cat = POT_NONE;
}

void fmd_pot_apply(fmd_t *md, unsigned atomkind1, unsigned atomkind2, fmd_pot_t *pot)
{
    // create the pottable if doesn't exist
    if (md->potsys.pottable == NULL) pottable_create(md);

    if (pot->cat == POT_EAM_ALLOY)
    {
        // find the local indices
        unsigned loc1, loc2;

        loc1 = _fmd_pot_eam_find_iloc(md, pot->data, atomkind1);
        loc2 = _fmd_pot_eam_find_iloc(md, pot->data, atomkind2);

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
    md->potsys.pottable[atomkind1][atomkind2].cat =
      md->potsys.pottable[atomkind2][atomkind1].cat = pot->cat;
    md->potsys.pottable[atomkind1][atomkind2].data =
      md->potsys.pottable[atomkind2][atomkind1].data = pot->data;

    if (md->potsys.potcats != NULL)
    {
        _fmd_list_free(md->potsys.potcats);
        md->potsys.potcats = NULL;
    }
}

static void process_potcats(fmd_t *md)
{
    md->potsys.eam_applied = false;

    // fill hybridpasses with zeros ( = flase )
    memset(md->potsys.hybridpasses, 0, sizeof(md->potsys.hybridpasses));

    list_t *potcats = md->potsys.potcats;
    while (potcats != NULL)
    {
        switch (*(potcat_t *)(potcats->data))
        {
            case POT_MORSE:
                md->potsys.hybridpasses[0] = true;
                break;

            case POT_LJ_6_12:
                md->potsys.hybridpasses[0] = true;
                break;

            case POT_EAM_ALLOY:
                md->potsys.hybridpasses[0] = true;
                md->potsys.hybridpasses[1] = true;
                md->potsys.eam_applied = true;
                break;
        }

        potcats = potcats->next;
    }
}

static int potcat_compare(const void *a, const void *b)
{
    if ( *( (potcat_t *)a ) == *( (potcat_t *)b ) )
        return 0;
    else
        return 1;
}

void _fmd_pot_update_and_process_potcats(fmd_t *md)
{
    // TO-DO: error should be handled here
    assert(md->potsys.pottable != NULL);

    if (md->potsys.potcats != NULL) return;  /* return if potcats is already updated */

    md->potsys.potcats_num = 0;

    for (unsigned i=0; i < md->potsys.atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
        {
            potpair_t *potpair = &md->potsys.pottable[i][j];

            // TO-DO: error should be handled here
            assert(potpair->cat != POT_NONE);

            // add the potcat to potcats list, if isn't already included there
            if (_fmd_list_find_custom(md->potsys.potcats, &potpair->cat, potcat_compare) == NULL)
            {
                md->potsys.potcats_num++;
                md->potsys.potcats = _fmd_list_prepend(md->potsys.potcats, &potpair->cat);
            }
        }

    process_potcats(md);
}

fmd_real_t _fmd_pot_get_largest_cutoff(potsys_t *ps)
{
    // TO-DO: error should be handled here
    assert(ps->pottable != NULL);

    fmd_real_t max = 0.0;

    for (unsigned i=0; i < ps->atomkinds_num; i++)
        for (unsigned j=0; j <= i; j++)
        {
            potpair_t *pair = &ps->pottable[i][j];
            fmd_real_t r;

            switch (pair->cat)
            {
                case POT_MORSE:
                    r = ((morse_t *)pair->data)->cutoff_sqr;
                    break;

                case POT_LJ_6_12:
                    r = ((LJ_6_12_t *)pair->data)->cutoff_sqr;
                    break;

                case POT_EAM_ALLOY:
                    r = ((eam_t *)pair->data)->cutoff_sqr;
                    break;

                default:
                    assert(0);  /* handle error */
            }

            if (r > max) max = r;
        }

    return sqrt(max);
}
