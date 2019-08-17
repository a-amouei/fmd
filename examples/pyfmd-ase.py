"""
   pyfmd-ase.py - an example showing how to use PyFMD in practice

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
   This example demonstrate how to build a structure using PyFMD interface for ASE and merge them together.
"""

from ase.visualize import view
import numpy as np
from ase.build import bulk
from ase.build.supercells import make_supercell

from pyfmd import *

# Create a super cell via ASE package
def make_supercell_cu():
    crys = bulk('Cu', 'fcc', a=3.6316, orthorhombic=True)
    P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) * 10
    return make_supercell(crys, P)

# Create molecular frames 1 (group-id = 0)
mf1 = MolecularFrame('Cu crystal 1').import_from(make_supercell_cu()).set_group_id(0)
# view(mf1.export())

# Create molecular frames 2 (group-id = 1)
mf2 = MolecularFrame('Cu crystal 2').import_from(make_supercell_cu()).set_group_id(1)
# view(mf2.export())

# Clone mf2 and translate atoms along given direction
mf2_moved = mf2.translate_atoms((30,-15,-15), inplace=False)

# Define region
def sphere(x, y, z):
    """sphere region"""
    if np.sqrt(x**2+y**2+z**2) < 10.0:
        return True
    else:
        return False

# Translate and select atoms within predefined region
mf1_sphere = mf1.translate_atoms((-10,-10,-10), inplace=False).select(Select().region(sphere))
# view(mf1_sphere.export())

# Merge two molecular frames
mf = mf1_sphere + mf2_moved
mf.translate_atoms((50, 50, 50))

# Define box and/or state header
box = Box((150, 100, 100))
state = StateHeader(simulation_time=0.0, number_of_atoms=mf.number_of_atoms, box=box, pbc=(False, False, False))
mf.set_sections([box, state])

# print(mf)
mf.write("system.stt")  # This state file can base used FMD simulation.
mf.read("system.stt")
view(mf.export())
