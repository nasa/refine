# -*- Makefile -*-

lib_LIBRARIES = \
	librefine.a \
	libFAKEGeomExtras.a \
	libFAKEGeomBox.a \
	libFAKEGeomPlate.a \
	libFAKEGeomTrap.a \
	libFAKEGeomDPWwbi.a \
	libFAKEGeomOM6v.a

librefine_a_CFLAGS = @SDKheader@
librefine_a_SOURCES = \
	refine_defs.h \
	adj.h adj.c \
	line.h line.c \
	queue.h queue.c \
	sort.h sort.c \
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
	gridfortran.h gridfortran.c

libFAKEGeomExtras_a_CFLAGS = @SDKheader@
libFAKEGeomExtras_a_SOURCES = FAKEGeom.h FAKEGeomExtras.c

libFAKEGeomBox_a_CFLAGS = @SDKheader@
libFAKEGeomBox_a_SOURCES = FAKEGeom.h FAKEGeomBox.c

libFAKEGeomPlate_a_CFLAGS = @SDKheader@
libFAKEGeomPlate_a_SOURCES = FAKEGeom.h FAKEGeomPlate.c

libFAKEGeomTrap_a_CFLAGS = @SDKheader@
libFAKEGeomTrap_a_SOURCES = FAKEGeom.h FAKEGeomTrap.c

libFAKEGeomDPWwbi_a_CFLAGS = @SDKheader@
libFAKEGeomDPWwbi_a_SOURCES = FAKEGeom.h FAKEGeomDPWwbi.c

libFAKEGeomOM6v_a_CFLAGS = @SDKheader@
libFAKEGeomOM6v_a_SOURCES = FAKEGeom.h FAKEGeomOM6v.c

if BUILD_EXECUTABLES

if BUILD_SDK_DEPENDENTS

lib_LIBRARIES += libMesherX.a

libMesherX_a_CFLAGS = @SDKheader@
libMesherX_a_SOURCES = \
	adj.h adj.c \
	near.h near.c \
	intersect.h intersect.c \
	line.h line.c \
	queue.h queue.c \
	sort.h sort.c \
	grid.h grid.c \
	gridmath.h gridmath.c \
	gridshape.h gridshape.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	gridmove.h gridmove.c \
	gridmpi.h gridmpi.c \
	gridfiller.h gridfiller.c \
	layer.h layer.c \
	mesherxInit.h mesherxInit.c \
	mesherxRebuild.h mesherxRebuild.c \
	mesherx.h mesherx.c

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

endif # BUILD_SDK_DEPENDENTS

EXTRA_DIST = mesherx_SyJet.c \
	makeRubyExtension.rb RubyExtensionBuilder.rb \
	adjMT.rb adj_ruby.c \
	nearMT.rb near_ruby.c \
	intersectMT.rb intersect_ruby.c \
	lineMT.rb line_ruby.c \
	sortMT.rb sort_ruby.c \
	queueMT.rb queue_ruby.c \
	gridMT.rb grid_ruby.c \
	gridmathMT.rb gridmath_ruby.c \
	gridmetricMT.rb gridmetric_ruby.c \
	gridshapeMT.rb gridshape_ruby.c \
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.h FAKEGeom.c \
	gridmpiMT.rb gridmpi_ruby.c \
	gridmoveMT.rb gridmove_ruby.c \
	layerMT.rb layer_ruby.c

TESTS = adjMT.rb nearMT.rb intersectMT.rb lineMT.rb sortMT.rb queueMT.rb gridMT.rb gridmathMT.rb gridmetricMT.rb gridshapeMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb gridmpiMT.rb gridmoveMT.rb layerMT.rb

bin_PROGRAMS =

if BUILD_SDK_DEPENDENTS

all_CFLAGS = @SDKheader@
refineLib  = libMesherX.a @SDKlibrary@ -lUGAPI -lGoolache -lCADGeom-CAPRI -lMeat
MesherXLib = libMesherX.a @SDKlibrary@ -lUGAPI -lGoolache -lCADGeom-CAPRI -lMeat
all_LDFLAGS = @dynamicFlag@
hla_sources = \
	MesherX_HLA.c \
	mesherx_hla.h mesherx_SyJet.c \
	layerNormalExtra.h layerNormalExtra.c \
	geometricStretch.h geometricStretch.c \
	SyJetSpacingControl.c

LastLib = @Xlibrary@ -lX11 -lm @ARCHlibrary@

if BUILD_CAPRI2

bin_PROGRAMS += refine MesherX MesherX_HLA
refine_CFLAGS = $(all_CFLAGS)
refine_SOURCES = refine_main.c
refine_LDADD = $(refineLib) @CAPRIlibrary@ $(LastLib)
refine_LDFLAGS = $(all_LDFLAGS)

MesherX_CFLAGS = $(all_CFLAGS)
MesherX_SOURCES = MesherX.c
MesherX_LDADD = $(MesherXLib) @CAPRIlibrary@ $(LastLib)
MesherX_LDFLAGS = $(all_LDFLAGS)

MesherX_HLA_CFLAGS = $(all_CFLAGS)
MesherX_HLA_SOURCES = $(hla_sources)
MesherX_HLA_LDADD = $(MesherXLib) @CAPRIlibrary@ $(LastLib)
MesherX_HLA_LDFLAGS = $(all_LDFLAGS)

bin_PROGRAMS += restartTest
restartTest_CFLAGS = $(all_CFLAGS)
restartTest_SOURCES = restartTest.c
restartTest_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	 @CAPRIlibrary@ $(LastLib)
restartTest_LDFLAGS = $(all_LDFLAGS)

endif # BUILD_CAPRI2

if BUILD_NATIVE_KERNEL

bin_PROGRAMS += refineNative MesherXNative MesherX_HLANative
refineNative_CFLAGS = $(all_CFLAGS)
refineNative_SOURCES = refine_main.c
refineNative_LDADD = $(refineLib) @CAPRINativelibrary@ $(LastLib)
refineNative_LDFLAGS = $(all_LDFLAGS)

MesherXNative_CFLAGS = $(all_CFLAGS)
MesherXNative_SOURCES = MesherX.c
MesherXNative_LDADD = $(MesherXLib) @CAPRINativelibrary@ $(LastLib)
MesherXNative_LDFLAGS = $(all_LDFLAGS)

MesherX_HLANative_CFLAGS = $(all_CFLAGS)
MesherX_HLANative_SOURCES = $(hla_sources)
MesherX_HLANative_LDADD = $(MesherXLib) @CAPRINativelibrary@ $(LastLib)
MesherX_HLANative_LDFLAGS = $(all_LDFLAGS)

bin_PROGRAMS += restartTestNative
restartTestNative_CFLAGS = $(all_CFLAGS)
restartTestNative_SOURCES = restartTest.c
restartTestNative_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	@CAPRINativelibrary@ $(LastLib)
restartTestNative_LDFLAGS = $(all_LDFLAGS)

bin_PROGRAMS += unitNative
unitNative_CFLAGS = $(all_CFLAGS)
unitNative_SOURCES = unit_main.c
unitNative_LDADD = $(refineLib) @CAPRINativelibrary@ $(LastLib)
unitNative_LDFLAGS = $(all_LDFLAGS)

if BUILD_PROE_KERNEL

bin_PROGRAMS += librefine-Linux.Pro.dll
librefine_Linux_Pro_dll_CFLAGS = $(all_CFLAGS) -DPROE_MAIN
librefine_Linux_Pro_dll_SOURCES = proemain.c
librefine_Linux_Pro_dll_LDADD = $(refineLib) @CAPRIProElibrary@ $(LastLib)
librefine_Linux_Pro_dll_LDFLAGS = $(all_LDFLAGS) -shared

bin_PROGRAMS += libMesherX-Linux.Pro.dll
libMesherX_Linux_Pro_dll_CFLAGS = $(all_CFLAGS) -DPROE_MAIN
libMesherX_Linux_Pro_dll_SOURCES = MesherX.c
libMesherX_Linux_Pro_dll_LDADD = $(MesherXLib) @CAPRIProElibrary@ $(LastLib)
libMesherX_Linux_Pro_dll_LDFLAGS = $(all_LDFLAGS) -shared

endif # BUILD_PROE_KERNEL

if BUILD_PARASOLID_KERNEL
bin_PROGRAMS += refineParasolid MesherXParasolid MesherX_HLAParasolid
refineParasolid_CFLAGS = $(all_CFLAGS)
refineParasolid_SOURCES = refine_main.c
refineParasolid_LDADD = $(refineLib) @CAPRIParasolidlibrary@ $(LastLib)
refineParasolid_LDFLAGS = $(all_LDFLAGS)

MesherXParasolid_CFLAGS = $(all_CFLAGS)
MesherXParasolid_SOURCES = MesherX.c
MesherXParasolid_LDADD = $(MesherXLib) @CAPRIParasolidlibrary@ $(LastLib)
MesherXParasolid_LDFLAGS = $(all_LDFLAGS)

MesherX_HLAParasolid_CFLAGS = $(all_CFLAGS)
MesherX_HLAParasolid_SOURCES = $(hla_sources)
MesherX_HLAParasolid_LDADD = $(MesherXLib) @CAPRIParasolidlibrary@ $(LastLib)
MesherX_HLAParasolid_LDFLAGS = $(all_LDFLAGS)

bin_PROGRAMS += restartTestParasolid
restartTestParasolid_CFLAGS = $(all_CFLAGS)
restartTestParasolid_SOURCES = restartTest.c
restartTestParasolid_LDADD = @SDKlibrary@ -lCADGeom-CAPRI -lMeat \
	@CAPRIParasolidlibrary@ $(LastLib)
restartTestParasolid_LDFLAGS = $(all_LDFLAGS)

bin_PROGRAMS += unitParasolid
unitParasolid_CFLAGS = $(all_CFLAGS)
unitParasolid_SOURCES = unit_main.c
unitParasolid_LDADD = $(refineLib) @CAPRIParasolidlibrary@ $(LastLib)
unitParasolid_LDFLAGS = $(all_LDFLAGS)

endif # BUILD_PARASOLID_KERNEL

if BUILD_OPENCASCADE_KERNEL
bin_PROGRAMS += refineOpenCASCADE MesherXOpenCASCADE MesherX_HLAOpenCASCADE
refineOpenCASCADE_CFLAGS = $(all_CFLAGS)
refineOpenCASCADE_SOURCES = refine_main.c
refineOpenCASCADE_LDADD = $(refineLib) @CAPRIOpenCASCADElibrary@ $(LastLib)
refineOpenCASCADE_LDFLAGS = $(all_LDFLAGS)

MesherXOpenCASCADE_CFLAGS = $(all_CFLAGS)
MesherXOpenCASCADE_SOURCES = MesherX.c
MesherXOpenCASCADE_LDADD = $(MesherXLib) @CAPRIOpenCASCADElibrary@ $(LastLib)
MesherXOpenCASCADE_LDFLAGS = $(all_LDFLAGS)

MesherX_HLAOpenCASCADE_CFLAGS = $(all_CFLAGS)
MesherX_HLAOpenCASCADE_SOURCES = $(hla_sources)
MesherX_HLAOpenCASCADE_LDADD = $(MesherXLib) @CAPRIOpenCASCADElibrary@ $(LastLib)
MesherX_HLAOpenCASCADE_LDFLAGS = $(all_LDFLAGS)
endif # BUILD_OPENCASCADE_KERNEL

endif # BUILD_NATIVE_KERNEL

endif # BUILD_SDK_DEPENDENTS

endif # BUILD_EXECUTABLES

