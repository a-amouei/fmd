/*
   03_ttm.c - an example showing how to use FMD in practice

   This program simulates interaction of short laser pulse with copper.
   It shows how to define a turi and perform a hybrid TTM-MD simulation
   with FMD.
*/

/* Assuming that FMD is already installed, this example can be compiled by

   $ gcc 03_ttm.c -lfmd -o 03_ttm.x

   and can be executed by

   $ mpirun -n 2 ./03_ttm.x
*/

#include <fmd.h>
#include <stdio.h>

fmd_handle_t turi, field_Te, field_Tl, field_n;
fmd_handle_t timer1;
FILE *fp_n, *fp_Te, *fp_Tl;

void saveRealArray(FILE *fp, fmd_utriple_t ut, fmd_array3_t a3)
{
    for (int i=0; i < ut[2]; i++)
        fprintf(fp, "%.2f ", ((fmd_real_t ***)a3)[0][0][i]);

    fprintf(fp, "\n");
    fflush(fp);
}

void saveUnsignedArray(FILE *fp, fmd_utriple_t ut, fmd_array3_t a3)
{
    for (int i=0; i < ut[2]; i++)
        fprintf(fp, "%u ", ((unsigned ***)a3)[0][0][i]);

    fprintf(fp, "\n");
    fflush(fp);
}

void handleEvents(fmd_t *md, fmd_event_t event, void *usp, fmd_params_t *params)
{
    switch (event)
    {
        case FMD_EVENT_TIMER_TICK: ;

            fmd_handle_t timer = ((fmd_event_params_timer_tick_t *)params)->timer;

            if (timer == timer1)
            {
                // report some quantities if the event is caused by timer1
                fmd_io_printf(md, "%f\t%f\n", fmd_dync_getTime(md),
                                              fmd_matt_getTotalEnergy(md));
            }

            break;

        case FMD_EVENT_FIELD_UPDATE: ;

            fmd_event_params_field_update_t *p = (fmd_event_params_field_update_t *)params;

            if ( p->turi == turi && (p->field == field_n || p->field == field_Te || p->field == field_Tl) )
            {
                fmd_utriple_t ut;
                fmd_array3_t a3;

                fmd_array3s_t *a3s = fmd_field_getArray(md, turi, p->field, &a3, ut);

                if (fmd_proc_isRoot(md))
                {
                    if (p->field == field_Te)
                        saveRealArray(fp_Te, ut, a3);
                    else if (p->field == field_Tl)
                        saveRealArray(fp_Tl, ut, a3);
                    else if (p->field == field_n)
                        saveUnsignedArray(fp_n, ut, a3);
                }

                fmd_array3s_free(a3s);
            }
    }
}

int main(int argc, char *argv[])
{
    fmd_t *md;
    double lp = 3.6316;       /* lattice parameter of copper */

    md = fmd_create();

    fmd_box_setSize(md, 10 * 10, 25 * lp, 40 * lp);

    fmd_box_setPBC(md, true, true, false);

    fmd_box_setSubdomains(md, 1, 1, 2);

    if (fmd_proc_isRoot(md))
    {
        fp_Te = fopen("Te.dat", "w");
        fp_Tl = fopen("Tl.dat", "w");
        fp_n = fopen("n.dat", "w");
    }

    // we have copper atoms
    fmd_string_t name[1] = {"Cu"};
    double mass[1] = {63.546};
    fmd_matt_setAtomKinds(md, 1, name, mass);

    // load the EAM file into memory; can be called only after fmd_box_setSubDomains()
    fmd_pot_t *pot = fmd_pot_eam_alloy_load(md, "../potentials/Cu01.eam.alloy");

    // apply the EAM potential
    fmd_pot_apply(md, 0, 0, pot);

    // make an fcc Cu cuboid at a given position and with a given size
    fmd_matt_makeCuboidFCC(md, 0.0, 0.0, 5.0*lp, 10, 10, 30, lp, 0, 0, 300.0);

    fmd_io_printf(md, "equilibrating the object...\n");
    fmd_dync_equilibrate(md, FMD_GROUP_ALL, .1, 2e-3, 2e-2, 300.0);

    fmd_setEventHandler(md, NULL, handleEvents);

    timer1 = fmd_timer_makeSimple(md, 0.0, 0.002, -1.0);

    turi = fmd_turi_add(md, FMD_TURI_TTM_TYPE1, 1, 1, 20, 0.0, -1.);
    field_Te = fmd_field_find(md, turi, FMD_FIELD_TTM_TE);
    field_Tl = fmd_field_find(md, turi, FMD_FIELD_TEMPERATURE);
    field_n = fmd_field_find(md, turi, FMD_FIELD_NUMBER);

    // set the parameters and quantities for TTM solver of type 1
    // the values are in SI unit system
    fmd_ttm_heat_capacity_linear_t C;
    C.gamma = 96.6;
    fmd_ttm_setHeatCapacity(md, turi, C);

    fmd_ttm_setHeatConductivity(md, turi, 400.);

    fmd_ttm_setCouplingFactor(md, turi, 1e17);

    fmd_ttm_setElectronTemperature(md, turi, 300.0);

    fmd_ttm_setTimestepRatio(md, turi, 200);

    fmd_ttm_setCellActivationFraction(md, turi, 0.1);

    fmd_ttm_laser_gaussian_t laser;
    laser.fluence = 4e4;
    laser.reflectance = 0.85;
    laser.duration = 100e-15;
    laser.t0 = 1e-12;
    laser.AbsorptionDepth = 14e-9;
    fmd_ttm_setLaserSource(md, turi, laser);

    fmd_io_printf(md, "start...\n");

    fmd_dync_integrate(md, FMD_GROUP_ALL, 2.0, 2e-3);

    if (fmd_proc_isRoot(md))
    {
        fclose(fp_Te);
        fclose(fp_Tl);
        fclose(fp_n);
    }

    fmd_free(md);

    return 0;
}
