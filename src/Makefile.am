# -*- Makefile -*-

# This file is included by the various $ARCH/Makefile.am

AUTOMAKE_OPTIONS = 1.7

lib_LIBRARIES = libMesherX.a librefine.a libOctree_Stub.a
libMesherX_a_CFLAGS = -I../src @SDKheader@
libMesherX_a_SOURCES = \
	adj.c \
	near.c \
	intersect.c \
	line.c \
	queue.c \
	sort.c \
	grid.c \
	gridmath.c \
	gridmetric.c \
	gridinsert.c \
	gridswap.c \
	gridcad.c \
	gridfiller.c \
	layer.c \
	mesherxInit.c \
	mesherxRebuild.c \
	mesherx.c

librefine_a_CFLAGS = -I../src @SDKheader@
librefine_a_SOURCES = \
	adj.c \
	near.c \
	line.c \
	queue.c \
	sort.c \
	grid.c \
	gridmath.c \
	gridmetric.c \
	gridinsert.c \
	gridswap.c \
	gridcad.c \
	gridmpi.c \
	gridfortran.c \
	FAKEGeomBox.c

libOctree_Stub_a_CFLAGS = -I../src @SDKheader@
libOctree_Stub_a_SOURCES = Octree_Stub.c

all_CFLAGS = -I../src @SDKheader@
refineLib = libMesherX.a @SDKlibrary@ -lCADGeom -lGoolache -lMeat
MesherXLib = libMesherX.a @SDKlibrary@ -lCADGeom -lGoolache libOctree_Stub.a @FELISAlibrary@ -ldl -lMeat
all_LDFLAGS = @dynamicFlag@
hla_sources = \
	MesherX_HLA.c \
	mesherx_hla.c \
	layerNormalExtra.c \
	geometricStretch.c \
	SyJetSpacingControl.c

mesherx_hla.c : mesherx_SyJet.c
	ln -s mesherx_SyJet.c mesherx_hla.c

LastLib = @Xlibrary@ -lX11 -lm @ARCHlibrary@

bin_PROGRAMS =

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
endif

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
endif

endif
