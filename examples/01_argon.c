/*
   01_argon.c - an example showing how to use FMD in practice

   Written in 2019 by the FMD authors

   To the extent possible under law, the author(s) have dedicated all
   copyright and related and neighboring rights to the current file to
   the public domain worldwide. This file is distributed without
   any warranty.

   You should have received a copy of the CC0 Public Domain Dedication
   along with this file. If not,
   see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

/* Assuming that FMD is already installed, this example can be compiled by

   $ gcc 01_argon.c -lfmd -lm -O3 -o 01_argon.x

   and executed by

   $ mpirun -n 2 ./01_argon.x
*/

#include <math.h>
#include <fmd.h>

int main(int argc, char *argv[])
{
    fmd_t *md;

    // create an FMD instance
    md = fmd_create();

    // set size of the simulation box (in Angstrom)
    double latticeParameter = 5.26;
    fmd_box_setSize(md, 10*latticeParameter, 10*latticeParameter, 10*latticeParameter);

    // set periodic boundary conditions in three dimensions (0 = no PBC)
    fmd_box_setPBC(md, 1, 1, 1);

    // partition the simulation box into subdomains for MPI-based parallel computation
    fmd_box_setSubDomains(md, 1, 2, 1);

    /* sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
       the function fmd_proc_isMD() can be called only after fmd_box_setSubDomains() */
    if (! fmd_proc_isMD(md))
    {
        fmd_free(md, 1);
        return 0;
    }

    // let's have only argon atoms
    fmd_string_t name[1] = {"Ar"};
    double mass[1] = {39.948};
    fmd_pot_setAtomKinds(md, 1, name, mass);

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
    fmd_matt_makeCuboidFCC(md, 0.0, 0.0, 0.0, 10, 10, 10, latticeParameter, 0, 0);

    // distribute the matter among subdomains
    fmd_matt_distribute(md);

    // set time step to 2 femtoseconds
    fmd_dync_setTimeStep(md, 2e-3);

    // let us simulate for 1.0 picoseconds
    double final_time = 1.0;

    // set where to save output files (default = current directory)
    //fmd_io_setSaveDirectory(md, "output/");

    // let configurations be saved as XYZ files
    fmd_io_setSaveConfigMode(md, FMD_SCM_XYZ_PARTICLESNUM);

    // set Berendsen thermostat parameter
    fmd_dync_setBerendsenThermostatParameter(md, 2e-2);

    // compute forces for the first time
    fmd_dync_updateForces(md);

    // the time loop starts here
    // fmd_dync_getTime() returns current internal time of the FMD instance
    while (fmd_dync_getTime(md) < final_time)
    {
        // save configuration every 40 femtoseconds
        if (fmod(fmd_dync_getTime(md), 0.04) < fmd_dync_getTimeStep(md))
            fmd_matt_saveConfiguration(md);

        // report some quantities every time step
        fmd_io_printf(md, "%f\t%f\t%e\n", fmd_dync_getTime(md),
                                          fmd_matt_getGlobalTemperature(md),
                                          fmd_matt_getTotalEnergy(md));

        // take first step of velocity Verlet integrator
        fmd_dync_velocityVerlet_takeFirstStep(md, 1);

        // compute forces
        fmd_dync_updateForces(md);

        // take last step of velocity Verlet integrator
        fmd_dync_velocityVerlet_takeLastStep(md);

        // increase internal time by one time step
        fmd_dync_incTime(md);
    }
    // end of the time loop

    // save system's final state in a file
    fmd_io_saveState(md, "state0.stt");

    // another report
    fmd_io_printf(md, "The run took about %.3f seconds to finish.\n", fmd_proc_getWallTime(md));

    // release memory taken for the FMD instance (including subdomain and all particles)
    // also finalize MPI
    fmd_free(md, 1);

    return 0;
}
