/*
   02_collision.c - an example showing how to use FMD in practice

   This program simulates a collision between two systems of atoms (Cu and Ar)
*/

/* Assuming that FMD is already installed, this example can be compiled by

   $ gcc 02_collision.c -lfmd -o 02_collision.x

   and executed by

   $ mpirun -n 2 ./02_collision.x
*/

#include <math.h>
#include <fmd.h>
#include <stdio.h>

typedef struct {fmd_handle_t timer1, timer2;} handles_t;

void handleEvents(fmd_t *md, fmd_event_t event, void *usp, fmd_params_t *params)
{
    switch (event)
    {
        case FMD_EVENT_TIMER_TICK: ;

            handles_t *handles = (handles_t *)usp;
            fmd_handle_t timer = ((fmd_event_params_timer_tick_t *)params)->timer;

            if (timer == handles->timer1)
            {
                // report some quantities if the event is caused by timer1
                fmd_io_printf(md, "%f\t%e\n", fmd_dync_getTime(md),
                                              fmd_matt_getTotalEnergy(md));

                fmd_rtriple_t p;

                fmd_matt_getMomentum(md, p);

                fmd_io_printf(md, "%f\t%f\t%f\n", p[0], p[1], p[2]);
            }
            else if (timer == handles->timer2)
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

    // create an FMD instance
    md = fmd_create();

    // set size of the simulation box (in Angstrom)
    double lx, ly, lz;
    fmd_box_setSize(md, lx=250.0, ly=250.0, lz=250.0);

    // partition the simulation box into subdomains for MPI-based parallel computation
    if (!fmd_box_setSubdomains(md, 1, 2, 1))
    { // sometimes the user launches more processes than the number of subdomains; they're not needed here!
        fmd_free(md);
        return 0;
    }

    // let's have copper and argon atoms
    fmd_string_t names[2] = {"Cu", "Ar"};
    double masses[2] = {63.546, 39.948};
    fmd_matt_setAtomKinds(md, 2, names, masses);

    // load the EAM file into memory; can be called only after fmd_box_setSubDomains()
    fmd_pot_t *pot = fmd_pot_eam_alloy_load(md, "../potentials/Cu01.eam.alloy");

    // apply the EAM potential for Cu
    fmd_pot_apply(md, 0, 0, pot);

    // use 12-6 Lennard-Jones potentials for Ar-Ar and Cu-Ar interactions
    fmd_pot_lj_apply(md, 1, 1, 3.40, 0.0104, 2.5*3.40);
    fmd_pot_lj_apply(md, 0, 1, 2.87, 0.0652, 2.5*2.87);

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
    fmd_matt_makeCuboidFCC(md, x0, y0, z0, cusize, cusize, cusize, lp0, 0, 0, 40.0);

    // add an fcc Ar cuboid with a different groupID
    fmd_matt_makeCuboidFCC(md, x1, y1, z1, cusize, cusize, cusize, lp1, 1, 1, 40.0);

    // set where to save output files (default = current directory)
    //fmd_io_setSaveDirectory(md, "output/");

    // save configurations as XYZ files
    fmd_io_setSaveConfigMode(md, FMD_SCM_XYZ_ATOMSNUM);

    // equilibrate the two colliding objects
    fmd_io_printf(md, "equilibrating the copper object...\n");
    fmd_dync_equilibrate(md, 0, 1.0, 2e-3, 2e-2, 40.0);
    fmd_io_printf(md, "equilibrating the argon object...\n");
    fmd_dync_equilibrate(md, 1, 1.0, 2e-3, 2e-2, 40.0);

    // add some center-of-mass velocity to the atoms of the objects (groups 0 and 1)
    fmd_matt_addVelocity(md, 0, +8., 0., 0.);
    fmd_matt_addVelocity(md, 1, -8., 0., 0.);

    handles_t handles;

    // assign an event handler to the FMD instance
    fmd_setEventHandler(md, &handles, handleEvents);

    // make two simple timers
    handles.timer1 = fmd_timer_makeSimple(md, 0.0, 0.05, -1.0);
    handles.timer2 = fmd_timer_makeSimple(md, 0.0, 0.06, -1.0);

    // simulate for 6.5 picoseconds, with timesteps of 2 fs
    fmd_dync_integrate(md, FMD_GROUP_ALL, 6.5, 2e-3);

    // save system's final state in a file
    //fmd_io_saveState(md, "state0.stt");

    fmd_io_printf(md, "The run took about %.3f seconds to finish.\n", fmd_proc_getWallTime(md));

    // release memory taken for the FMD instance (including subdomain and all particles)
    fmd_free(md);

    return 0;
}
