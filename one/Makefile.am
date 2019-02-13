
EXTRA_DIST = \
	makeRubyExtension.rb RubyExtensionBuilder.rb \
	adjMT.rb adj_ruby.c \
	lineMT.rb line_ruby.c \
	sortMT.rb sort_ruby.c \
	tableauMT.rb tableau_ruby.c \
	planMT.rb plan_ruby.c \
	queueMT.rb queue_ruby.c \
	gridMT.rb grid_ruby.c \
	gridmathMT.rb gridmath_ruby.c \
	gridmetricMT.rb gridmetric_ruby.c \
	gridshapeMT.rb gridshape_ruby.c \
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.c \
	gridmpiMT.rb gridmpi_ruby.c \
	gridmoveMT.rb gridmove_ruby.c \
	faux_input

lib_LIBRARIES =

lib_LIBRARIES += librefine.a

librefine_a_CFLAGS = @SDKheader@
librefine_a_SOURCES = \
	refine_defs.h \
	adj.h adj.c \
	line.h line.c \
	queue.h queue.c \
	sort.h sort.c \
	tableau.h tableau.c \
	plan.h plan.c \
	grid.h grid.c \
	gridmath.h gridmath.c \
	gridshape.h gridshape.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	gridgeom.h gridgeom.c \
	gridmpi.h gridmpi.c \
	gridmove.h gridmove.c \
	gridfortran.h gridfortran.c \
	unravel_fortran.h unravel_fortran.c

lib_LIBRARIES += libFAUXGeom.a

libFAUXGeom_a_CFLAGS = @SDKheader@
libFAUXGeom_a_SOURCES = FAKEGeom.h FAUXGeom.c FAKEGeomExtras.c

noinst_PROGRAMS =

noinst_PROGRAMS += refine-faux-test
refine_faux_test_SOURCES = FAUXTest.c
refine_faux_test_LDADD = libFAUXGeom.a librefine.a -lm

