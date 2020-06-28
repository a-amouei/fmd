/*
   03_timer_argon.c - an example showing how to use FMD in practice

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

   $ gcc 03_timer_argon.c -lfmd -O3 -o 03_timer_argon.x

   and can be executed by

   $ mpirun -n 2 ./03_timer_argon.x
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

    // create an fmd instance
    md = fmd_create();

    // assign an event handler to the instance
    fmd_setEventHandler(md, handleEvents);

    // make two simple timers (timer #0 and timer #1)
    fmd_timer_makeSimple(md, 0.0, 0.05, -1.0);
    fmd_timer_makeSimple(md, 0.0, 0.04, -1.0);

    // set size of the simulation box (in Angstrom)
    double latticeParameter = 5.26;
    fmd_box_setSize(md, 10*latticeParameter, 10*latticeParameter, 10*latticeParameter);

    // set periodic boundary conditions in three dimensions
    fmd_box_setPBC(md, FMD_TRUE, FMD_TRUE, FMD_TRUE);

    // partition the simulation box into subdomains for MPI-based parallel computation
    fmd_box_setSubDomains(md, 1, 2, 1);

    /* sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
       the function fmd_proc_isMD() can be called only after fmd_box_setSubDomains() */
    if (! fmd_proc_isMD(md))
    {
        fmd_free(md);
        return 0;
    }

    // let's have only argon atoms
    fmd_string_t name[1] = {"Ar"};
    double mass[1] = {39.948};
    fmd_matt_setAtomKinds(md, 1, name, mass);

    // use a 12-6 Lennard-Jones potential for Argon atoms
    double sigma = 3.4, epsilon = 0.0104;
    double cutoff = 2.5 * sigma;
    fmd_pot_lj_apply(md, 0, 0, sigma, epsilon, cutoff);

    // Here you can use a Morse potential for Argon atoms instead
    //double D0 = 0.010177, alpha = 1.253, r0 = 4.13, cutoff = 8.5;
    //fmd_pot_morse_apply(md, 0, 0, D0, alpha, r0, cutoff);

    // create the box grid
    fmd_box_createGrid(md, cutoff);

    // make an fcc cuboid at a given position and with a given size
    fmd_matt_makeCuboidFCC(md, 0.0, 0.0, 0.0, 10, 10, 10, latticeParameter, 0, 0);

    // distribute the matter among subdomains
    fmd_matt_distribute(md);

    // equilibrate for 1.0 picoseconds with a time step of 2 femtoseconds
    // to reach a temperature of 100 Kelvins
    fmd_dync_equilibrate(md, 0, 1.001, 2e-3, 2e-2, 100.0);

    // save system's final state in a file
    fmd_io_saveState(md, "state0.stt");

    // another report
    fmd_io_printf(md, "The run took about %.3f seconds to finish.\n", fmd_proc_getWallTime(md));

    // release memory taken for the fmd instance (including subdomain and all particles)
    fmd_free(md);

    return 0;
}
