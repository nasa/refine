bin_PROGRAMS = refine
refine_SOURCES = \
	master_header.h \
	adj.h adj.c \
	gridStruct.h \
	grid.h grid.c \
	gridswap.h gridswap.c \
	main.c
refine_LDADD = @CADGeom@ @Goolache@ @MeatLib@ @CAPRILib@
EXTRA_DIST = \
	adjMT.rb adj_ruby.c \
	gridMT.rb grid_ruby.c \
	gridswapMT.rb gridswap_ruby.c
TESTS = adjMT.rb gridMT.rb gridswapMT.rb
