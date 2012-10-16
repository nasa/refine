# -*- Makefile -*-

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

TESTS = adjMT.rb lineMT.rb sortMT.rb tableauMT.rb planMT.rb queueMT.rb gridMT.rb gridmathMT.rb gridmetricMT.rb gridshapeMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb gridmpiMT.rb gridmoveMT.rb

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

bin_PROGRAMS =

if BUILD_EXECUTABLES

bin_PROGRAMS += refine-convert
refine_convert_SOURCES = refine_convert.c
refine_convert_LDADD = librefine.a libFAUXGeom.a -lm

bin_PROGRAMS += refine-dsmc
refine_dsmc_SOURCES = refine_dsmc.c
refine_dsmc_LDADD = librefine.a libFAUXGeom.a -lm

bin_PROGRAMS += refine-px
refine_px_SOURCES = refine_px.c
refine_px_LDADD = librefine.a libFAUXGeom.a -lm

if BUILD_SDK_DEPENDENTS

all_cflags = @SDKheader@
refineLib  = librefine.a @SDKlibrary@ -lGoolache -lUGAPI -lCADGeom-CAPRI -lMeat
all_ldflags = @dynamicFlag@

LastLib = @Xlibrary@ -lX11 -lm @ARCHlibrary@

bin_PROGRAMS += refine
refine_CFLAGS = $(all_cflags)
refine_SOURCES = refine_main.c
refine_LDADD = $(refineLib) @CAPRIlibrary@ $(LastLib)
refine_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += restartTest
restartTest_CFLAGS = $(all_cflags)
restartTest_SOURCES = restartTest.c
restartTest_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	 @CAPRIlibrary@ $(LastLib)
restartTest_LDFLAGS = $(all_ldflags)

endif # BUILD_SDK_DEPENDENTS

endif # BUILD_EXECUTABLES
