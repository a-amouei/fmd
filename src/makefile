
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

PREFIX = /usr
CC = mpicc
SHAREDLIB = libfmd.so
CFLAGS = -fPIC -Wall -Wno-unused-result -O3 -fopenmp
LFLAGS = -fopenmp -lgsl -lgslcblas

OBJS ::= $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all shared clean install uninstall depend

all: shared

shared: $(SHAREDLIB)

$(SHAREDLIB): $(OBJS)
	$(CC) $(OBJS) -shared -o $@ $(LFLAGS)

include makefile.d

depend:
	gcc -MM *.c > makefile.d

clean:
	rm -f $(OBJS) $(SHAREDLIB)

install: $(SHAREDLIB)
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp fmd.h $(DESTDIR)$(PREFIX)/include/
	cp $(SHAREDLIB) $(DESTDIR)$(PREFIX)/lib/

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/include/fmd.h
	rm -f $(DESTDIR)$(PREFIX)/lib/$(SHAREDLIB)
