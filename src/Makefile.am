# -*- Makefile -*-

# This file is included by the various $ARCH/Makefile.am

AUTOMAKE_OPTIONS = 1.7

lib_LIBRARIES = \
	librefine.a \
	libFAKEGeomPlate.a \
	libFAKEGeomBox.a \
	libFAKEGeomTrap.a

librefine_a_CFLAGS = @SDKheader@
librefine_a_SOURCES = \
	master_header.h \
	adj.h adj.c \
	near.h near.c \
	line.h line.c \
	queue.h queue.c \
	sort.h sort.c \
	grid.h grid.c \
	gridmath.h gridmath.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	gridmpi.h gridmpi.c \
	gridfortran.h gridfortran.c

libFAKEGeomBox_a_CFLAGS = @SDKheader@
libFAKEGeomBox_a_SOURCES = FAKEGeomBox.h FAKEGeomBox.c

libFAKEGeomTrap_a_CFLAGS = @SDKheader@
libFAKEGeomTrap_a_SOURCES = FAKEGeomTrap.h FAKEGeomTrap.c

libFAKEGeomPlate_a_CFLAGS = @SDKheader@
libFAKEGeomPlate_a_SOURCES = FAKEGeomPlate.h FAKEGeomPlate.c

if BUILD_SDK_DEPENDENTS

lib_LIBRARIES += libMesherX.a libOctree_Stub.a

libMesherX_a_CFLAGS = @SDKheader@
libMesherX_a_SOURCES = \
	master_header.h \
	adj.h adj.c \
	near.h near.c \
	intersect.h intersect.c \
	line.h line.c \
	queue.h queue.c \
	sort.h sort.c \
	grid.h grid.c \
	gridmath.h gridmath.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	gridfiller.h gridfiller.c \
	layer.h layer.c \
	mesherxInit.h mesherxInit.c \
	mesherxRebuild.h mesherxRebuild.c \
	mesherx.h mesherx.c

libOctree_Stub_a_CFLAGS = @SDKheader@
libOctree_Stub_a_SOURCES = Octree_Stub.c

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
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.h FAKEGeom.c \
	gridmpiMT.rb gridmpi_ruby.c \
	layerMT.rb layer_ruby.c

TESTS = adjMT.rb nearMT.rb intersectMT.rb lineMT.rb sortMT.rb queueMT.rb gridMT.rb gridmathMT.rb gridmetricMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb gridmpiMT.rb layerMT.rb

bin_PROGRAMS =

if BUILD_SDK_DEPENDENTS
all_CFLAGS = @SDKheader@
refineLib = libMesherX.a @SDKlibrary@ -lCADGeom-CAPRI -lGoolache -lMeat
MesherXLib = libMesherX.a @SDKlibrary@ -lCADGeom-CAPRI -lGoolache libOctree_Stub.a -lMeat
all_LDFLAGS = @dynamicFlag@
hla_sources = \
	MesherX_HLA.c \
	mesherx_hla.h mesherx_hla.c \
	layerNormalExtra.h layerNormalExtra.c \
	geometricStretch.h geometricStretch.c \
	SyJetSpacingControl.c

mesherx_hla.c : mesherx_SyJet.c
	( cd $(srcdir) && ln -s mesherx_SyJet.c mesherx_hla.c )

LastLib = @Xlibrary@ -lX11 -lm @ARCHlibrary@

if BUILD_NATIVE_KERNEL
bin_PROGRAMS += refineNative MesherXNative MesherX_HLANative
refineNative_CFLAGS = $(all_CFLAGS)
refineNative_SOURCES = main.c
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

if BUILD_PARASOLID_KERNEL
bin_PROGRAMS += refineParasolid MesherXParasolid MesherX_HLAParasolid
refineParasolid_CFLAGS = $(all_CFLAGS)
refineParasolid_SOURCES = main.c
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
endif # BUILD_PARASOLID_KERNEL

if BUILD_OPENCASCADE_KERNEL
bin_PROGRAMS += refineOpenCASCADE MesherXOpenCASCADE MesherX_HLAOpenCASCADE
refineOpenCASCADE_CFLAGS = $(all_CFLAGS)
refineOpenCASCADE_SOURCES = main.c
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