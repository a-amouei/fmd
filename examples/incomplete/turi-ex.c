/*
    - an example showing how to use FMD in practice

   Written in 2020 by the authors of FMD

   To the extent possible under law, the author(s) have dedicated all
   copyright and related and neighboring rights to the current file to
   the public domain worldwide. This file is distributed without
   any warranty.

   You should have received a copy of the CC0 Public Domain Dedication
   along with this file. If not,
   see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

/* Assuming that FMD is already installed, this example can be compiled by

   $ gcc turi-ex.c -lfmd -lm -O3 -o turi-ex.x

   and executed by

   $ mpirun -n 2 ./turi-ex.x
*/

#include <math.h>
#include <fmd.h>
#include <stdio.h>

fmd_handle_t turi, field;

void handleEvents(fmd_t *md, fmd_event_t event, fmd_event_params_t *params)
{
    static int counter = 0;

    switch (event)
    {
        case FMD_EVENT_FIELD_UPDATE: ;
            fmd_event_params_field_update_t *p = (fmd_event_params_field_update_t *)params;
            if (p->field == field && p->turi == turi)
            {
                char str[32];

                fmd_io_printf(md, "%f\t%s\n", fmd_dync_getTime(md), "number-density field updated!");

                sprintf(str, "out-%05d.h5", counter++);

                fmd_field_save_as_hdf5(md, turi, field, str);
            }
            else
                fmd_io_printf(md, "%f\t%s\n", fmd_dync_getTime(md), "another field updated!");

            break;
    }
}

int main(int argc, char *argv[])
{
    fmd_t *md;

    // create an FMD instance
    md = fmd_create();

    // assign an event handler to the instance
    fmd_setEventHandler(md, handleEvents);

    // set size of the simulation box (in Angstrom)
    double lx, ly, lz;
    fmd_box_setSize(md, lx=250.0, ly=250.0, lz=250.0);

    // set periodic boundary conditions in three dimensions
    fmd_box_setPBC(md, FMD_FALSE, FMD_FALSE, FMD_FALSE);

    // partition the simulation box into subdomains for MPI-based parallel computation
    fmd_box_setSubdomains(md, 1, 2, 1);

    /* sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
       the function fmd_proc_isMD() can be called only after fmd_box_setSubDomains() */
    if (! fmd_proc_isMD(md))
    {
        fmd_free(md);
        return 0;
    }

    // let's have copper and argon atoms
    fmd_string_t names[2] = {"Cu", "Ar"};
    double masses[2] = {63.546, 39.948};
    fmd_matt_setAtomKinds(md, 2, names, masses);

    // load the EAM file into memory; can be called only after fmd_box_setSubDomains()
    fmd_pot_t *pot = fmd_pot_eam_alloy_load(md, "../../potentials/Cu01.eam.alloy");

    // apply the EAM potential for Cu
    fmd_pot_apply(md, 0, 0, pot);

    // use 12-6 Lennard-Jones potentials for Ar-Ar and Cu-Ar interactions
    fmd_pot_lj_apply(md, 1, 1, 3.40, 0.0104, 2.5*3.40);
    fmd_pot_lj_apply(md, 0, 1, 2.87, 0.0652, 2.5*2.87);

    // create the grid
    fmd_box_createGrid(md, 2.5*3.40);

    // prepare some parameters
    double dc = 30.0;     // the initial distance between the colliding objects
    int cusize = 7;       // cusize x cusize x cusize = number of unit cells in each object
    double lp0 = 3.6316;  // lattice parameter of copper
    double lp1 = 5.26;    // lattice parameter of argon
    double x0 = (lx - dc) / 2 - cusize * lp0;
    double y0 = (ly - cusize * lp0) / 2;
    double z0 = (lz - cusize * lp0) / 2;
    double x1 = (lx + dc) / 2;
    double y1 = (ly - cusize * lp1) / 2;
    double z1 = (lz - cusize * lp1) / 2;

    // make an fcc Cu cuboid at a given position and with a given size
    fmd_matt_makeCuboidFCC(md, x0, y0, z0, cusize, cusize, cusize, lp0, 0, 0);

    // add an fcc Ar cuboid with a different groupID
    fmd_matt_makeCuboidFCC(md, x1, y1, z1, cusize, cusize, cusize, lp1, 1, 1);

    // distribute the matter among subdomains
    fmd_matt_distribute(md);

    // set time step to 2 femtoseconds
    fmd_dync_setTimeStep(md, 2e-3);

    // set where to save output files (default = current directory)
    //fmd_io_setSaveDirectory(md, "output/");

    // save configurations as XYZ files
    fmd_io_setSaveConfigMode(md, FMD_SCM_XYZ_PARTICLESNUM);

    // equilibrate the two colliding objects
    fmd_io_printf(md, "equilibrating the copper object...\n");
    fmd_dync_equilibrate(md, 0, 1., 2e-3, 2e-2, 40.0);
    fmd_io_printf(md, "equilibrating the argon object...\n");
    fmd_dync_equilibrate(md, 1, 1., 2e-3, 2e-2, 40.0);

    // add some center-of-mass velocity to the atoms of the objects (groups 0 and 1)
    fmd_matt_addVelocity(md, 0, +8., 0., 0.);
    fmd_matt_addVelocity(md, 1, -8., 0., 0.);

    turi = fmd_turi_add(md, FMD_TURI_CUSTOM, 10, 10, 10, 0.0, -1.);
    field = fmd_field_add(md, turi, FMD_FIELD_NUMBER_DENSITY, 0.1);

    // activate all groups for dynamics; -1 as a groupID means all groups
    fmd_matt_setActiveGroup(md, -1);

    // simulate for 6.5 picoseconds
    double final_time = 1.0;

    // compute forces for the first time
    fmd_dync_updateForces(md);

    // the time loop starts here
    // fmd_dync_getTime() returns current internal time of the FMD instance
    while (fmd_dync_getTime(md) < final_time)
    {
        // save configuration every 60 femtoseconds
        if (fmod(fmd_dync_getTime(md), 0.06) < fmd_dync_getTimeStep(md))
            fmd_matt_saveConfiguration(md);

        // report some quantities every time step
        /*fmd_io_printf(md, "%f\t%e\n", fmd_dync_getTime(md),
                                      fmd_matt_getTotalEnergy(md));*/

        // use velocity Verlet integrator: start step
        fmd_dync_VelocityVerlet_startStep(md, FMD_FALSE);

        // compute forces
        fmd_dync_updateForces(md);

        // use velocity Verlet integrator: finish step
        fmd_dync_VelocityVerlet_finishStep(md);

        // increase internal time by one time step
        fmd_dync_incTime(md);
    }
    // end of the time loop

    // save system's final state in a file
    //fmd_io_saveState(md, "state0.stt");

    // another report
    fmd_io_printf(md, "The run took about %.3f seconds to finish.\n", fmd_proc_getWallTime(md));

    // release memory taken for the FMD instance (including subdomain and all particles)
    fmd_free(md);

    return 0;
}
