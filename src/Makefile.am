lib_LIBRARIES = libMesherX.a
libMesherX_a_SOURCES = \
	master_header.h \
	adj.h adj.c \
	near.h near.c \
	intersect.h intersect.c \
	grid.h grid.c \
	gridmath.h gridmath.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	gridfiller.h gridfiller.c \
	layer.h layerStruct.h layer.c \
	mesherxInit.h mesherxInit.c \
	mesherxRebuild.h mesherxRebuild.c \
	mesher.h mesherx.c 

include_HEADERS = \
	master_header.h \
	mesherx.h \
	adj.h \
	near.h \
	intersect.h \
	grid.h \
	layer.h

bin_PROGRAMS = refine MesherX MesherX_HLA
refine_SOURCES = main.c
refine_LDADD = libMesherX.a @CADGeom@ @Goolache@ @MeatLib@ @CAPRILib@ -lm

MesherX_SOURCES = MesherX.c
MesherX_LDADD = libMesherX.a @CADGeom@ @Goolache@ @FELISALib@ @MeatLib@ @CAPRILib@ -lm
MesherX_LDFLAGS = -rdynamic

MesherX_HLA_SOURCES = MesherX_HLA.c \
	mesherx_hla.c mesherx_hla.h \
	layerNormalExtra.h layerNormalExtra.c \
	geometricStretch.h geometricStretch.c
MesherX_HLA_LDADD = libMesherX.a @CADGeom@ @Goolache@ @FELISALib@ @MeatLib@ @CAPRILib@ -lm
MesherX_HLA_LDFLAGS = -rdynamic

EXTRA_DIST = makeRubyExtension.rb RubyExtensionBuilder.rb \
	adjMT.rb adj_ruby.c \
	nearMT.rb near_ruby.c \
	intersectMT.rb intersect_ruby.c \
	gridMT.rb grid_ruby.c \
	gridmathMT.rb gridmath_ruby.c \
	gridmetricMT.rb gridmetric_ruby.c \
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.h FAKEGeom.c \
	layerMT.rb layer_ruby.c

TESTS = adjMT.rb nearMT.rb intersectMT.rb gridMT.rb gridmetricMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb layerMT.rb MesherX
