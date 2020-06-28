array.o: array.c array.h config.h
base.o: base.c base.h config.h potential.h types.h md_ghost.h forces.h \
 timer.h molecule.h
cspline.o: cspline.c cspline.h config.h
eam.o: eam.c eam.h config.h forces.h types.h cspline.h base.h potential.h \
 list.h
forces.o: forces.c forces.h config.h eam.h types.h cspline.h lj.h morse.h \
 base.h potential.h md_ghost.h list.h molecule.h
list.o: list.c list.h config.h
lj.o: lj.c lj.h config.h potential.h types.h base.h list.h forces.h
md_ghost.o: md_ghost.c base.h config.h potential.h types.h md_ghost.h
molecule.o: molecule.c molecule.h config.h types.h potential.h base.h \
 array.h list.h
morse.o: morse.c morse.h config.h base.h potential.h types.h forces.h \
 list.h
potential.o: potential.c potential.h config.h types.h base.h array.h \
 list.h eam.h forces.h cspline.h molecule.h
structure.o: structure.c base.h config.h potential.h types.h molecule.h \
 array.h
timer.o: timer.c timer.h config.h types.h base.h potential.h
