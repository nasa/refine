bin_PROGRAMS = refine
refine_SOURCES = \
	master_header.h \
	adj.h adj.c \
	gridStruct.h \
	grid.h grid.c \
	gridmetric.h gridmetric.c \
	gridinsert.h gridinsert.c \
	gridswap.h gridswap.c \
	gridcad.h gridcad.c \
	main.c
refine_LDADD = @CADGeom@ @Goolache@ @MeatLib@ @CAPRILib@
EXTRA_DIST = makeRubyExtension.rb \
	adjMT.rb adj_ruby.c \
	gridMT.rb grid_ruby.c \
	gridmetricMT.rb gridmetric_ruby.c \
	gridinsertMT.rb gridinsert_ruby.c \
	gridswapMT.rb gridswap_ruby.c \
	gridcadMT.rb gridcad_ruby.c \
	FAKEGeom.h FAKEGeom.c
TESTS = adjMT.rb gridMT.rb gridmetricMT.rb gridinsertMT.rb gridswapMT.rb gridcadMT.rb
