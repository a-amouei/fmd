
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
# (rarely needed if the makefile can find mpicc)
MPI_HEADER =

# the path to the "lib(64)" directory of MPI (the directory having "libmpi.so")
# (rarely needed if the makefile can find mpicc)
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

PREFIX ?= /usr/local
BUILDPATH = ../build
LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include
PKGCONFIGDIR = $(LIBDIR)/pkgconfig
VERSION = 0.4.0

# set CC

ifeq ($(strip $(CC)),)
	CC = mpicc
endif

# set MPI_LIB, I_MPI, L_MPI

ifeq ($(strip $(MPI_HEADER)),)
	I_MPI =
else
	I_MPI = -I$(MPI_HEADER)
endif

ifeq ($(strip $(MPI_LIB)),)
	MPICC_SHOW := $(shell $(CC) -show 2>/dev/null)
	MPI_LIB := $(shell echo $(MPICC_SHOW) | sed -n 's/.*-L\([^ ]*\).*/\1/p')
endif

ifeq ($(strip $(MPI_LIB)),)
	L_MPI =
else
	L_MPI = -L$(MPI_LIB)
endif

# set HDF5_LIB, I_HDF5, L_HDF5

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

# set GSL_LIB, I_GSL, L_GSL

ifeq ($(strip $(GSL_HEADER)),)
	GSL_CFLAGS := $(shell pkg-config --cflags gsl 2>/dev/null)
	GSL_HEADER := $(shell echo $(GSL_CFLAGS) | sed -n 's/.*-I\([^ ]*\).*/\1/p')
endif

ifeq ($(strip $(GSL_HEADER)),)
	I_GSL =
else
	I_GSL = -I$(GSL_HEADER)
endif

ifeq ($(strip $(GSL_LIB)),)
	GSL_LIBS := $(shell pkg-config --libs gsl 2>/dev/null)
	GSL_LIB := $(shell echo $(GSL_LIBS) | sed -n 's/.*-L\([^ ]*\).*/\1/p')
endif

ifeq ($(strip $(GSL_LIB)),)
	L_GSL =
else
	L_GSL = -L$(GSL_LIB)
endif

# Variables for the .pc file
NAME = fmd
DESCRIPTION = Free Molecular Dynamics
PKGCPATH = $(PKGCONFIGDIR)/$(NAME).pc
# End of variables for the .pc file

HDF5_RPATH = $(HDF5_LIB)
MPI_RPATH = $(MPI_LIB)
GSL_RPATH = $(GSL_LIB)

LIBNAME = libfmd
SHAREDLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a

CFLAGS = -fPIC -std=c11 -pedantic -O3 -fopenmp $(I_HDF5) $(I_MPI) $(I_GSL)
LFLAGS = -fopenmp -Wl,-rpath,$(HDF5_RPATH):$(MPI_RPATH):$(GSL_RPATH) $(L_MPI) $(L_HDF5) $(L_GSL) -lgsl -lgslcblas -lhdf5

OBJS ::= $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all shared clean install uninstall depend

all: shared static

shared: $(SHAREDLIB)

$(SHAREDLIB): $(OBJS)
	$(CC) -shared -o $@ $(OBJS) $(LFLAGS)

static: $(STATICLIB)

$(STATICLIB): $(OBJS)
	ar rcs $(STATICLIB) $(OBJS)

include makefile.d

depend:
	gcc -MM *.c > makefile.d

clean:
	@rm -f *.o $(SHAREDLIB) $(STATICLIB)
	@rm -f $(BUILDPATH)/$(SHAREDLIB)
	@rm -f $(BUILDPATH)/$(STATICLIB)
	@rm -f $(BUILDPATH)/fmd.h
	@rm -fd $(BUILDPATH)
	@echo "Cleaned!"

install: $(SHAREDLIB) $(STATICLIB)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(INCLUDEDIR)
	@cp fmd.h $(INCLUDEDIR)
	@cp $(SHAREDLIB) $(LIBDIR)
	@cp $(STATICLIB) $(LIBDIR)
	@mkdir -p $(PKGCONFIGDIR)
	@echo 'prefix=$(PREFIX)' > $(PKGCPATH)
	@echo 'exec_prefix=$${prefix}' >> $(PKGCPATH)
	@echo 'libdir=${LIBDIR}' >> $(PKGCPATH)
	@echo 'includedir=${INCLUDEDIR}\n' >> $(PKGCPATH)
	@echo 'Name: $(NAME)' >> $(PKGCPATH)
	@echo 'Description: $(DESCRIPTION)' >> $(PKGCPATH)
	@echo 'Version: $(VERSION)' >> $(PKGCPATH)
	@echo 'Requires.private: mpi gsl hdf5' >> $(PKGCPATH)
	@echo 'Libs: -L$${libdir} -lfmd' >> $(PKGCPATH)
	@echo 'Libs.private: -fopenmp' >> $(PKGCPATH)
	@echo 'Cflags: -I$${includedir}' >> $(PKGCPATH)
	@ldconfig
	@echo "The library was installed successfully!"

build: $(SHAREDLIB)
	@mkdir -p $(BUILDPATH)
	@cp fmd.h $(BUILDPATH)
	@cp $(SHAREDLIB) $(BUILDPATH)
	@cp $(STATICLIB) $(BUILDPATH)
	@echo "The library was built successfully!"

uninstall:
	@rm -f $(INCLUDEDIR)/fmd.h
	@rm -f $(LIBDIR)/$(SHAREDLIB)
	@rm -f $(LIBDIR)/$(STATICLIB)
	@rm -f $(PKGCPATH)
	@ldconfig
	@echo "The library was uninstalled successfully!"