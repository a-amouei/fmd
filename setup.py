#   setup.py: This file is part of Free Molecular Dynamics
#
#   Copyright (C) 2019 Hossein Ghorbanfekr
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys
from setuptools import setup

if sys.version_info[:2] < (3, 5):
    print("Python >= 3.5 is required.")
    sys.exit(-1)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pyfmd',
#    version='0.1',
    description='PyFMD provides an object-oriented interface for interacting with the core part of FMD.',
    license='GNU GPL v3',
    url='https://github.com/a-amouei/fmd',
    keywords="fmd molecular-dynamics physics",
    platforms=['Unix'],
    long_description=long_description,
    long_description_content_type='text/markdown',
#    author='',
#    author_email='',
    packages=['pyfmd'],  # 'src'
    install_requires=['numpy', 'ase', 'periodictable']  #
)
