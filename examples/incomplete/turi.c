/*
   turi.c - an example showing how to use FMD in practice

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

   $ gcc turi.c -lfmd -O3 -o turi.x

   and can be executed by

   $ mpirun -n 2 ./turi.x
*/

#include <fmd.h>

int main(int argc, char *argv[])
{
    fmd_t *md;

    md = fmd_create();

    fmd_box_setSize(md, 50., 50., 50.);

    fmd_box_setSubDomains(md, 2, 2, 2);

    fmd_box_createGrid(md, 1.0);

    fmd_turi_add(md, FMD_TURI_CUSTOM, 5, 5, 5);

    fmd_free(md);
    return 0;
}
