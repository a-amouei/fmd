array.o: array.c array.h config.h types.h general.h
cell.o: cell.c cell.h config.h types.h fmd-private.h subdomain.h array.h \
 potential.h events.h h5.h misc.h general.h
cspline.o: cspline.c cspline.h config.h types.h general.h
eam.o: eam.c eam.h config.h forces.h types.h cspline.h fmd-private.h \
 subdomain.h array.h potential.h cell.h events.h h5.h misc.h list.h \
 general.h
forces.o: forces.c forces.h config.h fmd-private.h types.h subdomain.h \
 array.h potential.h cell.h events.h h5.h eam.h cspline.h lj.h morse.h \
 misc.h md-ghost.h list.h molecule.h general.h turi.h ttm.h
general.o: general.c general.h config.h types.h
h5.o: h5.c h5.h config.h types.h fmd-private.h subdomain.h array.h \
 potential.h cell.h events.h turi.h general.h misc.h
integrators.o: integrators.c fmd-private.h config.h types.h subdomain.h \
 array.h potential.h cell.h events.h h5.h matter.h misc.h turi.h forces.h \
 timer.h general.h ttm.h
list.o: list.c list.h config.h general.h types.h
lj.o: lj.c lj.h config.h types.h fmd-private.h subdomain.h array.h \
 potential.h cell.h events.h h5.h misc.h list.h forces.h general.h
matter.o: matter.c matter.h config.h types.h fmd-private.h subdomain.h \
 array.h potential.h cell.h events.h h5.h general.h molecule.h
md-ghost.o: md-ghost.c fmd-private.h config.h types.h subdomain.h array.h \
 potential.h cell.h events.h h5.h misc.h md-ghost.h general.h
misc.o: misc.c matter.h config.h types.h misc.h fmd-private.h subdomain.h \
 array.h potential.h cell.h events.h h5.h md-ghost.h timer.h turi.h \
 general.h
molecule.o: molecule.c molecule.h config.h types.h fmd-private.h \
 subdomain.h array.h potential.h cell.h events.h h5.h misc.h list.h \
 general.h
morse.o: morse.c morse.h config.h types.h fmd-private.h subdomain.h \
 array.h potential.h cell.h events.h h5.h misc.h forces.h list.h \
 general.h
potential.o: potential.c potential.h config.h types.h fmd-private.h \
 subdomain.h array.h cell.h events.h h5.h misc.h list.h eam.h forces.h \
 cspline.h molecule.h general.h morse.h lj.h
structure.o: structure.c fmd-private.h config.h types.h subdomain.h \
 array.h potential.h cell.h events.h h5.h misc.h molecule.h general.h
subdomain.o: subdomain.c subdomain.h config.h types.h array.h \
 fmd-private.h potential.h cell.h events.h h5.h misc.h general.h
timer.o: timer.c timer.h config.h types.h fmd-private.h subdomain.h \
 array.h potential.h cell.h events.h h5.h misc.h general.h
ttm.o: ttm.c ttm.h config.h types.h array.h fmd-private.h subdomain.h \
 potential.h cell.h events.h h5.h matter.h general.h turi.h turi-ghost.h
turi-ghost.o: turi-ghost.c turi-ghost.h config.h types.h fmd-private.h \
 subdomain.h array.h potential.h cell.h events.h h5.h misc.h general.h \
 turi.h
turi.o: turi.c turi.h config.h types.h array.h fmd-private.h subdomain.h \
 potential.h cell.h events.h h5.h misc.h timer.h general.h ttm.h
version.o: version.c types.h config.h general.h
