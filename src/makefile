
#  makefile: This file is part of Free Molecular Dynamics

#  Copyright (C) 2019 Arham Amouye Foumani

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#################################################################################
# If needed, set the following paths according to your system configuration. If
# one or more are left empty, this makefile, the compiler, and the linker try to
# find the correct paths. However, if this fails and you receive error messages,
# or if there are multiple versions of one or more of the following libraries
# on disk, you may consider setting some/all of the paths explicitly.
#################################################################################

# the path to the file "mpicc" (not the directory having it)
CC =

# the path to the "include" directory of MPI
MPI_HEADER =

# the path to the "lib(64)" directory of MPI (the directory having "libmpi.so")
MPI_LIB =

# the path to the "include" directory of HDF5
HDF5_HEADER =

# the path to the "lib(64)" directory of HDF5 (the directory having "libhdf5.so")
HDF5_LIB =

# the path to the "include" directory of GSL
GSL_HEADER =

# the path to the "lib(64)" directory of GSL (the directory having "libgsl.so")
GSL_LIB =
#################################################################################

INSTALL_PREFIX = /usr
BUILDPATH = ../build
SHAREDLIB = libfmd.so

ifeq ($(strip $(CC)),)
	CC = mpicc
endif

ifeq ($(strip $(MPI_HEADER)),)
	I_MPI =
else
	I_MPI = -I$(MPI_HEADER)
endif

ifeq ($(strip $(MPI_LIB)),)
	L_MPI =
else
	L_MPI = -L$(MPI_HEADER)
endif

ifeq ($(strip $(HDF5_HEADER)),)
	HDF5_HEADER := $(shell command -v h5cc >/dev/null 2>&1 && h5cc -show | sed -n 's/.*-I\([^ ]*\).*/\1/p')
endif

ifeq ($(strip $(HDF5_HEADER)),)
	I_HDF5 =
else
	I_HDF5 = -I$(HDF5_HEADER)
endif

ifeq ($(strip $(HDF5_LIB)),)
	HDF5_LIB := $(shell command -v h5cc >/dev/null 2>&1 && h5cc -show | sed -n 's/.*-L\([^ ]*\).*/\1/p')
endif

ifeq ($(strip $(HDF5_LIB)),)
	L_HDF5 =
else
	L_HDF5 = -L$(HDF5_LIB)
endif

ifeq ($(strip $(GSL_HEADER)),)
	I_GSL =
else
	I_GSL = -I$(GSL_HEADER)
endif

ifeq ($(strip $(GSL_LIB)),)
	L_GSL =
else
	L_GSL = -L$(GSL_LIB)
endif

HDF5_RPATH = $(HDF5_LIB)
MPI_RPATH = $(MPI_LIB)
GSL_RPATH = $(GSL_LIB)

CFLAGS = -fPIC -std=c11 -O3 -fopenmp $(I_HDF5) $(I_MPI) $(I_GSL)
LFLAGS = -fopenmp -Wl,-rpath,$(HDF5_RPATH):$(MPI_RPATH):$(GSL_RPATH) $(L_MPI) $(L_HDF5) $(L_GSL) -lgsl -lgslcblas -lhdf5

OBJS ::= $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all shared clean install uninstall depend

all: shared

shared: $(SHAREDLIB)

$(SHAREDLIB): $(OBJS)
	$(CC) -shared -o $@ $(OBJS) $(LFLAGS)

include makefile.d

depend:
	gcc -MM *.c > makefile.d

clean:
	rm -f *.o $(SHAREDLIB)
	rm -f $(BUILDPATH)/$(SHAREDLIB)
	rm -f $(BUILDPATH)/fmd.h
	rm -fd $(BUILDPATH)

install: $(SHAREDLIB)
	mkdir -p $(DESTDIR)$(INSTALL_PREFIX)/lib
	mkdir -p $(DESTDIR)$(INSTALL_PREFIX)/include
	cp fmd.h $(DESTDIR)$(INSTALL_PREFIX)/include/
	cp $(SHAREDLIB) $(DESTDIR)$(INSTALL_PREFIX)/lib/
	ldconfig

build: $(SHAREDLIB)
	mkdir -p $(BUILDPATH)
	cp fmd.h $(BUILDPATH)
	cp $(SHAREDLIB) $(BUILDPATH)

uninstall:
	rm -f $(DESTDIR)$(INSTALL_PREFIX)/include/fmd.h
	rm -f $(DESTDIR)$(INSTALL_PREFIX)/lib/$(SHAREDLIB)
	ldconfig
