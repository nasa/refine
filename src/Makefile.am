lib_LIBRARIES = libMesherX.a
libMesherX_a_SOURCES = \
	master_header.h \
	adj.h adj.c \
	grid.h grid.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	gridfiller.h gridfiller.c \
	near.h near.c \
	layer.h layerStruct.h layer.c \
	mesherxInit.h mesherxInit.c \
	mesherxRebuild.h mesherxRebuild.c \
	mesher.h mesherx.c 

include_HEADERS = \
	master_header.h \
	mesherx.h \
	adj.h \
	grid.h \
	near.h \
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

EXTRA_DIST = makeRubyExtension.rb \
	adjMT.rb adj_ruby.c \
	gridMT.rb grid_ruby.c \
	gridmetricMT.rb gridmetric_ruby.c \
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.h FAKEGeom.c \
	nearMT.rb near_ruby.c \
	layerMT.rb layer_ruby.c

TESTS = adjMT.rb gridMT.rb gridmetricMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb layerMT.rb nearMT.rb MesherX
