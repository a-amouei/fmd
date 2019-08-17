"""
   02_collision.py - an example showing how to use PyFMD in practice

   Written in 2019 by the FMD authors

   To the extent possible under law, the author(s) have dedicated all
   copyright and related and neighboring rights to the current file to
   the public domain worldwide. This file is distributed without
   any warranty.

   You should have received a copy of the CC0 Public Domain Dedication
   along with this file. If not,
   see <http://creativecommons.org/publicdomain/zero/1.0/>.
"""

"""
   This example simulates collision between two nanocrystals (Ar & Cu) in NVE ensemble.
   Before running it, make sure that PyFMD and the core part of FMD are both installed.

   How to run:

       $ mpirun -n 2 python 02_collision.py
"""

import sys
from pyfmd import Calculator

# create an fmd-system instance
md = Calculator()

# set size of the simulation box (in Angstrom)
lx, ly, lz = 250.0, 250.0, 250.0
md.box_size = lx, ly, lz

# set periodic boundary conditions in three dimensions (False = no PBC)
md.box_pbc = False, False, False

# partition the simulation box into subdomains for MPI-based parallel computation
md.subdomains = (1, 1, 2)

# sometimes the user launches more processes than the chosen number of subdomains; they're not needed here!
# 'is_process_md' can be called only after setting subdomains.
if not md.is_process_md:
    md.free_system(finalize_mpi=True)
    sys.exit()

# add copper and argon atoms
md.set_potential_atom_kinds({'Cu': 63.546, 'Ar': 39.948})  # symbols and masses (a.m.u)

# load the EAM file into memory; can be called only after setting subdomains
pot = md.load_potential_eam_alloy("../potentials/Cu01.eam.alloy")

# apply the EAM potential for Cu
md.apply_potential(0, 0, pot)

# use 12-6 Lennard-Jones potentials for Ar-Ar and Cu-Ar interactions
md.apply_potential_lj(1, 1, 3.40, 0.0104, 2.5*3.40)
md.apply_potential_lj(0, 1, 2.87, 0.0652, 2.5*2.87)

# create the box grid
md.create_box_grid(2.5*3.40)

# set the desired temperature (in Kelvin)
md.desired_temperature = 40.0

# prepare some parameters
dc = 30.0      # the initial distance between the colliding objects
cusize = 7     # cusize x cusize x cusize = number of unit cells in each object
lp0 = 3.6316   # lattice parameter of copper
lp1 = 5.26     # lattice parameter of argon
x0 = (lx - dc) / 2 - cusize * lp0
y0 = (ly - cusize * lp0) / 2
z0 = (lz - cusize * lp0) / 2
x1 = (lx + dc) / 2
y1 = (ly - cusize * lp1) / 2
z1 = (lz - cusize * lp1) / 2

# make an fcc Cu cuboid at a given position and with a given size
md.make_cuboid_fcc((x0, y0, z0), (cusize, cusize, cusize), lp0, 0, 0)

# add an fcc Ar cuboid with a different groupID
md.make_cuboid_fcc((x1, y1, z1), (cusize, cusize, cusize), lp1, 1, 1)

# distribute the matter among subdomains
md.material_distribute()

# set time step to 2 femtoseconds
md.time_step = 0.002

# set where to save output files (default = current directory)
md.io_directory = "../"

# let configurations be saved as XYZ files
md.set_io_config_mode(0)

# set Berendsen thermostat parameter
md.set_berendsen_thermostat_parameter(0.02)

# equilibrate the two clusters
print("equilibrating the copper object...\n")
md.equilibrate(0, 1.0, 2e-2)
print("equilibrating the argon object...\n")
md.equilibrate(1, 1.0, 2e-2)

# add some center-of-mass velocity to atoms of groups 0 and 1
md.add_velocity(0, (+8., 0., 0.))
md.add_velocity(1, (-8., 0., 0.))

# activate all groups for dynamics; -1 as a groupID means all groups
md.set_active_group(-1)

# compute forces for the first time
md.update_force()

final_step = 6.5
while md.time < final_step:

    # save configuration every 60 femtoseconds
    if md.get_time() % 0.06 < md.get_time_step():
        md.save_io_config()

    # report some quantities every time step
    # print (md.time, md.temperature, md.total_energy)

    # time integrator
    md.velocity_verlet(use_thermostat=False)  # NVE

    # increase internal time by one time step
    md.next_time_step()

# save system's final state in a file
md.save_io_state("state0.stt")

# another report
print("The run took about %.3f seconds to finish." % md.process_wall_time)

# release memory taken for potential and fmd-system
# del md
