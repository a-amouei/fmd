"""
   pyfmd-molframe.py - an example showing how to use PyFMD in practice

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
   This example demonstrate how to build a structure using PyFMD interface.
"""


from ase.visualize import view
import numpy as np

from pyfmd import MolecularFrame, Atom, Atoms, Box, Select

# Make and distribute atoms
atoms = Atoms()

L = [_ for _ in range(0, 20, 2)]
L_max = max(L)

for x in L:
    for y in L:
        for z in L:
            atom = Atom(symbol='X', position=(x, y, z), group_id=int(x + y + z) % 3)
            atoms.append(atom)

# Make a molecular frame
mf = MolecularFrame("Sample MF", atoms)

# print (mf)
print ("number of atoms (all):", mf.number_of_atoms)
# view(mf.export())  # draw atoms

# Select atoms with specific group-id
sel1 = Select().group_id([0, 1])
# view(mf.select(sel1).export())  # draw atoms

# Translate atoms and box
mf.translate_atoms((-L_max / 2, -L_max / 2, -L_max / 2))


# Define a region (sphere)
def sphere(x, y, z):
    if np.sqrt(x ** 2 + y ** 2 + z ** 2) < 8:
        return True
    else:
        return False


# Select atoms with a sphere region
sel2 = Select().region(sphere)
mf2 = mf.select(sel2)
# view(mf2.export())  # draw atoms

# logical operator between selections
# view(mf.select(sel2.AND(sel1)).export())  # draw atoms

# Make box section and add to molecular frame
box = Box((L_max, L_max, L_max))
mf2.set_sections(box)

# Set back atoms into center
mf2.translate_atoms((L_max / 2, L_max / 2, L_max / 2))

print ("number of atoms (sphere):", mf2.number_of_atoms)
view(mf2.export())  # draw atoms

# Save to file
mf2.write("sphere.xyz")
# mf2.write("sphere.stt")
