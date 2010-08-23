# -*- Makefile -*-

EXTRA_DIST = \
	makeRubyExtension.rb RubyExtensionBuilder.rb \
	adjMT.rb adj_ruby.c \
	nearMT.rb near_ruby.c \
	intersectMT.rb intersect_ruby.c \
	lineMT.rb line_ruby.c \
	sortMT.rb sort_ruby.c \
	tableauMT.rb tableau_ruby.c \
	planMT.rb plan_ruby.c \
	ringMT.rb ring_ruby.c \
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
	layerMT.rb layer_ruby.c \
	faux_input

TESTS = adjMT.rb nearMT.rb intersectMT.rb lineMT.rb sortMT.rb tableauMT.rb planMT.rb ringMT.rb queueMT.rb gridMT.rb gridmathMT.rb gridmetricMT.rb gridshapeMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb gridmpiMT.rb gridmoveMT.rb layerMT.rb

lib_LIBRARIES =

lib_LIBRARIES += librefine.a

librefine_a_CFLAGS = @SDKheader@
librefine_a_SOURCES = \
	refine_defs.h \
	adj.h adj.c \
	intersect.h intersect.c \
	line.h line.c \
	queue.h queue.c \
	sort.h sort.c \
	tableau.h tableau.c \
	plan.h plan.c \
	ring.h ring.c \
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
refine_px_LDADD = libMesherX.a libFAUXGeom.a -lm

if BUILD_SDK_DEPENDENTS

lib_LIBRARIES += libMesherX.a

mesherx_sources = \
	gridfiller.h gridfiller.c \
	layer.h layer.c \
	mesherxInit.h mesherxInit.c \
	mesherxRebuild.h mesherxRebuild.c \
	mesherx.h mesherx.c

libMesherX_a_CFLAGS = @SDKheader@
libMesherX_a_SOURCES = $(mesherx_sources)
libMesherX_a_include_HEADERS = \
	refine_defs.h \
	mesherx.h \
	layer.h \
	grid.h \
	near.h \
	adj.h \
	line.h \
	queue.h

libMesherX_a_includedir = @prefix@/include

all_cflags = @SDKheader@
refineLib  = libMesherX.a @SDKlibrary@ -lGoolache -lUGAPI -lCADGeom-CAPRI -lMeat
MesherXLib = libMesherX.a @SDKlibrary@ -lGoolache -lUGAPI -lCADGeom-CAPRI -lMeat
all_ldflags = @dynamicFlag@

LastLib = @Xlibrary@ -lX11 -lm @ARCHlibrary@

if BUILD_CAPRI2

bin_PROGRAMS += refine
refine_CFLAGS = $(all_cflags)
refine_SOURCES = refine_main.c
refine_LDADD = $(refineLib) @CAPRIlibrary@ $(LastLib)
refine_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += MesherX
MesherX_CFLAGS = $(all_cflags)
MesherX_SOURCES = mesherx_main.c
MesherX_LDADD = $(refineLib) $(MesherXLib) @CAPRIlibrary@ $(LastLib)
MesherX_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += restartTest
restartTest_CFLAGS = $(all_cflags)
restartTest_SOURCES = restartTest.c
restartTest_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	 @CAPRIlibrary@ $(LastLib)
restartTest_LDFLAGS = $(all_ldflags)

endif # BUILD_CAPRI2

if BUILD_NATIVE_KERNEL

bin_PROGRAMS += refineNative
refineNative_CFLAGS = $(all_cflags)
refineNative_SOURCES = refine_main.c
refineNative_LDADD = $(refineLib) @CAPRINativelibrary@ $(LastLib)
refineNative_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += MesherXNative
MesherXNative_CFLAGS = $(all_cflags)
MesherXNative_SOURCES = mesherx_main.c
MesherXNative_LDADD = $(refineLib) $(MesherXLib) @CAPRINativelibrary@ $(LastLib)
MesherXNative_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += restartTestNative
restartTestNative_CFLAGS = $(all_cflags)
restartTestNative_SOURCES = restartTest.c
restartTestNative_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	@CAPRINativelibrary@ $(LastLib)
restartTestNative_LDFLAGS = $(all_ldflags)

if BUILD_PROE_KERNEL

bin_PROGRAMS += librefine-Linux.Pro.dll
librefine_Linux_Pro_dll_CFLAGS = $(all_cflags) -DPROE_MAIN
librefine_Linux_Pro_dll_SOURCES = refine_main.c
librefine_Linux_Pro_dll_LDADD = $(refineLib) @CAPRIProElibrary@ $(LastLib)
librefine_Linux_Pro_dll_LDFLAGS = $(all_ldflags) -shared

bin_PROGRAMS += libMesherX-Linux.Pro.dll
libMesherX_Linux_Pro_dll_CFLAGS = $(all_cflags) -DPROE_MAIN
libMesherX_Linux_Pro_dll_SOURCES = mesherx_main.c
libMesherX_Linux_Pro_dll_LDADD = $(refineLib) $(MesherXLib) @CAPRIProElibrary@ $(LastLib)
libMesherX_Linux_Pro_dll_LDFLAGS = $(all_ldflags) -shared

endif # BUILD_PROE_KERNEL

if BUILD_PARASOLID_KERNEL
bin_PROGRAMS += refineParasolid
refineParasolid_CFLAGS = $(all_cflags)
refineParasolid_SOURCES = refine_main.c
refineParasolid_LDADD = $(refineLib) @CAPRIParasolidlibrary@ $(LastLib)
refineParasolid_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += MesherXParasolid
MesherXParasolid_CFLAGS = $(all_cflags)
MesherXParasolid_SOURCES = mesherx_main.c
MesherXParasolid_LDADD = $(refineLib) $(MesherXLib) @CAPRIParasolidlibrary@ $(LastLib)
MesherXParasolid_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += restartTestParasolid
restartTestParasolid_CFLAGS = $(all_cflags)
restartTestParasolid_SOURCES = restartTest.c
restartTestParasolid_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	@CAPRIParasolidlibrary@ $(LastLib)
restartTestParasolid_LDFLAGS = $(all_ldflags)

endif # BUILD_PARASOLID_KERNEL

if BUILD_OPENCASCADE_KERNEL
bin_PROGRAMS += refineOpenCASCADE
refineOpenCASCADE_CFLAGS = $(all_cflags)
refineOpenCASCADE_SOURCES = refine_main.c
refineOpenCASCADE_LDADD = $(refineLib) @CAPRIOpenCASCADElibrary@ $(LastLib)
refineOpenCASCADE_LDFLAGS = $(all_ldflags)

bin_PROGRAMS += MesherXOpenCASCADE
MesherXOpenCASCADE_CFLAGS = $(all_cflags)
MesherXOpenCASCADE_SOURCES = mesherx_main.c
MesherXOpenCASCADE_LDADD = $(refineLib) $(MesherXLib) @CAPRIOpenCASCADElibrary@ $(LastLib)
MesherXOpenCASCADE_LDFLAGS = $(all_ldflags)

endif # BUILD_OPENCASCADE_KERNEL

endif # BUILD_NATIVE_KERNEL

endif # BUILD_SDK_DEPENDENTS

endif # BUILD_EXECUTABLES

