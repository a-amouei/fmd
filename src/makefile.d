array.o: array.c array.h config.h types.h general.h
base.o: base.c base.h config.h potential.h types.h subdomain.h array.h \
 events.h h5.h md_ghost.h forces.h timer.h molecule.h turi.h general.h
cspline.o: cspline.c cspline.h config.h types.h
eam.o: eam.c eam.h config.h forces.h types.h cspline.h base.h potential.h \
 subdomain.h array.h events.h h5.h list.h general.h
forces.o: forces.c forces.h config.h eam.h types.h cspline.h lj.h morse.h \
 base.h potential.h subdomain.h array.h events.h h5.h md_ghost.h list.h \
 molecule.h general.h
general.o: general.c general.h config.h types.h
h5.o: h5.c h5.h config.h types.h array.h turi.h general.h base.h \
 potential.h subdomain.h events.h
list.o: list.c list.h config.h
lj.o: lj.c lj.h config.h types.h potential.h base.h subdomain.h array.h \
 events.h h5.h list.h forces.h general.h
md_ghost.o: md_ghost.c base.h config.h potential.h types.h subdomain.h \
 array.h events.h h5.h md_ghost.h general.h
molecule.o: molecule.c molecule.h config.h types.h potential.h base.h \
 subdomain.h array.h events.h h5.h list.h general.h
morse.o: morse.c morse.h config.h types.h base.h potential.h subdomain.h \
 array.h events.h h5.h forces.h list.h general.h
potential.o: potential.c potential.h config.h types.h base.h subdomain.h \
 array.h events.h h5.h list.h eam.h forces.h cspline.h molecule.h
structure.o: structure.c base.h config.h potential.h types.h subdomain.h \
 array.h events.h h5.h molecule.h general.h
subdomain.o: subdomain.c subdomain.h config.h types.h array.h base.h \
 potential.h events.h h5.h general.h
timer.o: timer.c timer.h config.h types.h base.h potential.h subdomain.h \
 array.h events.h h5.h
turi.o: turi.c turi.h config.h types.h array.h base.h potential.h \
 subdomain.h events.h h5.h timer.h general.h
