lib_LIBRARIES = libMesherX.a
libMesherX_a_SOURCES = \
	master_header.h \
	adj.h adj.c \
	gridStruct.h \
	grid.h grid.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	layer.h layer.c

include_HEADERS = \
	master_header.h \
	grid.h \
	layer.h

bin_PROGRAMS = refine
refine_SOURCES = main.c
refine_LDADD = libMesherX.a @CADGeom@ @Goolache@ @MeatLib@ @CAPRILib@

EXTRA_DIST = makeRubyExtension.rb \
	adjMT.rb adj_ruby.c \
	gridMT.rb grid_ruby.c \
	gridmetricMT.rb gridmetric_ruby.c \
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.h FAKEGeom.c \
	layerMT.rb layer_ruby.c

TESTS = adjMT.rb gridMT.rb gridmetricMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb layerMT.rb
