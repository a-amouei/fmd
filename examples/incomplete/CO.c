/*
   CO.c - an example showing how to use FMD in practice

   Written in 2019 by the authors of FMD

   To the extent possible under law, the author(s) have dedicated all
   copyright and related and neighboring rights to the current file to
   the public domain worldwide. This file is distributed without
   any warranty.

   You should have received a copy of the CC0 Public Domain Dedication
   along with this file. If not,
   see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

/* Assuming that FMD is already installed, this example can be compiled by

   $ gcc CO.c -lfmd -O3 -o CO.x

   and can be executed by

   $ mpirun -n 2 ./CO.x
*/

#include <fmd.h>

void handleEvents(fmd_t *md, fmd_event_t event, unsigned param)
{
    switch (event)
    {
        case FMD_EVENT_TIMERTICK:
            if (param == 0)
            {
                // report some quantities if the event is caused by timer #0
                fmd_io_printf(md, "%f\t%f\t%e\n", fmd_dync_getTime(md),
                                                  fmd_matt_getGlobalTemperature(md),
                                                  fmd_matt_getTotalEnergy(md));
            }
            else if (param == 1)
            {
                // save configuration if the event is caused by timer #1
                fmd_matt_saveConfiguration(md);
            }
            break;
    }
}

int main(int argc, char *argv[])
{
    fmd_t *md;

    md = fmd_create();

    // assign an event handler to the instance
    fmd_setEventHandler(md, handleEvents);

    // make two simple timers (timer #0 and timer #1)
    fmd_timer_makeSimple(md, 0.0, 0.05, -1.0);
    fmd_timer_makeSimple(md, 0.0, 0.04, -1.0);

    fmd_box_setSize(md, 50., 50., 50.);

    fmd_box_setPBC(md, 1, 1, 1);

    fmd_box_setSubDomains(md, 1, 1, 1);

    fmd_string_t atom_names[2] = {"C", "O"};
    double atom_masses[2] = {12.011, 15.999};
    fmd_matt_setAtomKinds(md, 2, atom_names, atom_masses);

    fmd_pot_lj_apply(md, 0, 0, 3.6, 2.7e-3, 2.5*3.6);
    fmd_pot_lj_apply(md, 1, 1, 3.0, 4.5e-3, 2.5*3.0);
    fmd_pot_lj_apply(md, 0, 1, 3.3, 3.5e-3, 2.5*3.3);

    double bond_coeffs[] = {58.0, 1.13};     // k (eV/Å^2) and r0 (Å)
    unsigned CO_bond = fmd_bond_addKind(md, FMD_BOND_HARMONIC, bond_coeffs);

    unsigned CO_atoms_kinds[2] = {0, 1};
    double CO_atoms_positions[2][3] = {{0., 0., 0.}, {1.13, 0., 0.}};
    unsigned CO_mol = fmd_molecule_addKind(md, "Carbon Monoxide", 2,
                                           CO_atoms_kinds, CO_atoms_positions);

    fmd_bond_apply(md, CO_bond, CO_mol, 0, 1);

    fmd_box_createGrid(md, 2.5*3.6);

    fmd_matt_scatterMolecule(md, CO_mol, 0., 0., 0., 50., 50., 50., 35, 0);

    fmd_matt_distribute(md);

    fmd_dync_equilibrate(md, 0, 100.001, 0.1e-3, 0.1e-2, 300.0);

    fmd_free(md);
    return 0;
}
