# -*- Makefile -*-

lib_LIBRARIES = libref2.a

libref2_a_SOURCES = \
	ref_adj.c ref_adj.h \
	ref_axi.c ref_axi.h \
	ref_cell.c ref_cell.h \
	ref_dict.c ref_dict.h \
	ref_defs.h \
	ref_edge.c ref_edge.h \
	ref_endian.h \
	ref_export.c ref_export.h \
	ref_face.c ref_face.h \
	ref_fixture.c ref_fixtusre.h \
	ref_fortran.c ref_fortran.h \
	ref_gather.c ref_gather.h \
	ref_grid.c ref_grid.h \
	ref_import.c ref_import.h \
	ref_list.h ref_list.c \
	ref_math.h ref_math.c \
	ref_malloc.h \
	ref_metric.h ref_metric.c \
	ref_migrate.h ref_migrate.c \
	ref_mpi.h ref_mpi.c \
	ref_node.h ref_node.c \
	ref_part.h ref_part.c \
	ref_quality.h ref_quality.c \
	ref_shard.c ref_shard.h \
	ref_sort.h ref_sort.c \
	ref_subdiv.h ref_subdiv.c \
	ref_swap.h ref_swap.c \
	ref_validation.h ref_validation.c
libref2_a_CFLAGS = @zoltan_INCLUDE@

bin_PROGRAMS =

bin_PROGRAMS += dpw5_faceid_prism
dpw5_faceid_prism_SOURCES = dpw5_faceid_prism.c
dpw5_faceid_prism_LDADD = libref2.a -lm

bin_PROGRAMS += dpw5_pole_prism
dpw5_pole_prism_SOURCES = dpw5_pole_prism.c
dpw5_pole_prism_LDADD = libref2.a -lm

bin_PROGRAMS += dpw5_prism_tet
dpw5_prism_tet_SOURCES = dpw5_prism_tet.c
dpw5_prism_tet_LDADD = libref2.a -lm

bin_PROGRAMS += projectable
projectable_SOURCES = projectable.c
projectable_LDADD = libref2.a -lm

bin_PROGRAMS += convert
convert_SOURCES = convert.c
convert_LDADD = libref2.a -lm

noinst_PROGRAMS =
TESTS =

TESTS += ref_adj_test
noinst_PROGRAMS += ref_adj_test
ref_adj_test_SOURCES = ref_adj_test.c
ref_adj_test_LDADD = libref2.a -lm

TESTS += ref_axi_test
noinst_PROGRAMS += ref_axi_test
ref_axi_test_SOURCES = ref_axi_test.c
ref_axi_test_LDADD = libref2.a -lm

TESTS += ref_cell_test
noinst_PROGRAMS += ref_cell_test
ref_cell_test_SOURCES = ref_cell_test.c
ref_cell_test_LDADD = libref2.a -lm

TESTS += ref_dict_test
noinst_PROGRAMS += ref_dict_test
ref_dict_test_SOURCES = ref_dict_test.c
ref_dict_test_LDADD = libref2.a -lm

TESTS += ref_edge_test
noinst_PROGRAMS += ref_edge_test
ref_edge_test_SOURCES = ref_edge_test.c
ref_edge_test_LDADD = libref2.a -lm

TESTS += ref_export_test
noinst_PROGRAMS += ref_export_test
ref_export_test_SOURCES = ref_export_test.c
ref_export_test_LDADD = libref2.a -lm

TESTS += ref_face_test
noinst_PROGRAMS += ref_face_test
ref_face_test_SOURCES = ref_face_test.c
ref_face_test_LDADD = libref2.a -lm

TESTS += ref_fixture_test
noinst_PROGRAMS += ref_fixture_test
ref_fixture_test_SOURCES = ref_fixture_test.c
ref_fixture_test_LDADD = libref2.a @zoltan_LDADD@ -lm

TESTS += ref_fortran_test
noinst_PROGRAMS += ref_fortran_test
ref_fortran_test_SOURCES = ref_fortran_test.c
ref_fortran_test_LDADD = libref2.a -lm

TESTS += ref_gather_test
noinst_PROGRAMS += ref_gather_test
ref_gather_test_SOURCES = ref_gather_test.c
ref_gather_test_LDADD = libref2.a @zoltan_LDADD@ -lm

TESTS += ref_grid_test
noinst_PROGRAMS += ref_grid_test
ref_grid_test_SOURCES = ref_grid_test.c
ref_grid_test_LDADD = libref2.a -lm

TESTS += ref_import_test
noinst_PROGRAMS += ref_import_test
ref_import_test_SOURCES = ref_import_test.c
ref_import_test_LDADD = libref2.a -lm

TESTS += ref_list_test
noinst_PROGRAMS += ref_list_test
ref_list_test_SOURCES = ref_list_test.c
ref_list_test_LDADD = libref2.a -lm

TESTS += ref_math_test
noinst_PROGRAMS += ref_math_test
ref_math_test_SOURCES = ref_math_test.c
ref_math_test_LDADD = libref2.a -lm

TESTS += ref_metric_test
noinst_PROGRAMS += ref_metric_test
ref_metric_test_SOURCES = ref_metric_test.c
ref_metric_test_LDADD = libref2.a -lm

TESTS += ref_migrate_test
noinst_PROGRAMS += ref_migrate_test
ref_migrate_test_SOURCES = ref_migrate_test.c
ref_migrate_test_LDADD = libref2.a @zoltan_LDADD@ -lm

TESTS += ref_mpi_test
noinst_PROGRAMS += ref_mpi_test
ref_mpi_test_SOURCES = ref_mpi_test.c
ref_mpi_test_LDADD = libref2.a -lm

TESTS += ref_node_test
noinst_PROGRAMS += ref_node_test
ref_node_test_SOURCES = ref_node_test.c
ref_node_test_LDADD = libref2.a -lm

TESTS += ref_part_test
noinst_PROGRAMS += ref_part_test
ref_part_test_SOURCES = ref_part_test.c
ref_part_test_LDADD = libref2.a @zoltan_LDADD@ -lm

TESTS += ref_quality_test
noinst_PROGRAMS += ref_quality_test
ref_quality_test_SOURCES = ref_quality_test.c
ref_quality_test_LDADD = libref2.a -lm

TESTS += ref_shard_test
noinst_PROGRAMS += ref_shard_test
ref_shard_test_SOURCES = ref_shard_test.c
ref_shard_test_LDADD = libref2.a -lm

TESTS += ref_sort_test
noinst_PROGRAMS += ref_sort_test
ref_sort_test_SOURCES = ref_sort_test.c
ref_sort_test_LDADD = libref2.a -lm

TESTS += ref_subdiv_test
noinst_PROGRAMS += ref_subdiv_test
ref_subdiv_test_SOURCES = ref_subdiv_test.c
ref_subdiv_test_LDADD = libref2.a -lm

TESTS += ref_swap_test
noinst_PROGRAMS += ref_swap_test
ref_swap_test_SOURCES = ref_swap_test.c
ref_swap_test_LDADD = libref2.a -lm

TESTS += ref_validation_test
noinst_PROGRAMS += ref_validation_test
ref_validation_test_SOURCES = ref_validation_test.c
ref_validation_test_LDADD = libref2.a -lm
