/*
   01_argon.c - an example showing how to use FMD in practice

   This program equilibrates a system of atoms.
*/

/* Assuming that FMD is already installed, this example can be compiled by

   $ gcc 01_argon.c -lfmd -O3 -o 01_argon.x

   and can be executed by

   $ mpirun -n 2 ./01_argon.x
*/

#include <fmd.h>

fmd_handle_t timer1, timer2;

void handleEvents(fmd_t *md, fmd_event_t event, fmd_params_t *params)
{
    switch (event)
    {
        case FMD_EVENT_TIMER_TICK: ;

            fmd_handle_t timer = ((fmd_event_params_timer_tick_t *)params)->timer;

            if (timer == timer1)
            {
                // report some quantities if the event is caused by timer1
                fmd_io_printf(md, "%f\t%f\t%e\n", fmd_dync_getTime(md),
                                                  fmd_matt_getTemperature(md),
                                                  fmd_matt_getTotalEnergy(md));
            }
            else if (timer == timer2)
            {
                // save configuration if the event is caused by timer2
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

    // make two simple timers
    timer1 = fmd_timer_makeSimple(md, 0.0, 0.05, -1.0);
    timer2 = fmd_timer_makeSimple(md, 0.0, 0.04, -1.0);

    // set size of the simulation box (in Angstrom)
    double LatticeParameter = 5.26;
    fmd_box_setSize(md, 10*LatticeParameter, 10*LatticeParameter, 10*LatticeParameter);

    // set periodic boundary conditions in three dimensions
    fmd_box_setPBC(md, true, true, true);

    // partition the simulation box into subdomains for MPI-based parallel computation
    fmd_box_setSubdomains(md, 1, 2, 1);

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

    // set the desired temperature (in Kelvin)
    fmd_matt_setDesiredTemperature(md, 100.0);

    // make an fcc cuboid at a given position and with a given size
    fmd_matt_makeCuboidFCC(md, 0.0, 0.0, 0.0, 10, 10, 10, LatticeParameter, 0, 0);

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
