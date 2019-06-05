
EXTRA_DIST = test.sh

lib_LIBRARIES = libref2.a

libref2_a_SOURCES = \
	ref_adapt.c \
	ref_adj.c \
	ref_agents.c \
	ref_args.c \
	ref_axi.c \
	ref_cavity.c \
	ref_cell.c \
	ref_cloud.c \
	ref_clump.c \
	ref_collapse.c \
	ref_comprow.c \
	ref_dict.c \
	ref_edge.c \
	ref_elast.c \
	ref_endian.h \
	ref_export.c \
	ref_face.c \
	ref_fixture.c \
	ref_fortran.c \
	ref_gather.c \
	ref_geom.c \
	ref_grid.c \
	ref_histogram.c \
	ref_html.c \
	ref_import.c \
	ref_inflate.c \
	ref_interp.c \
	ref_layer.c \
	ref_list.c \
	ref_math.c \
	ref_matrix.c \
	ref_metric.c \
	ref_migrate.c \
	ref_mpi.c \
	ref_node.c \
	ref_part.c \
	ref_phys.c \
	ref_recon.c \
	ref_search.c \
	ref_shard.c \
	ref_smooth.c \
	ref_sort.c \
	ref_split.c \
	ref_subdiv.c \
	ref_swap.c \
	ref_twod.c \
	ref_validation.c

include_HEADERS = ref_adapt.h ref_adj.h ref_agents.h ref_args.h ref_axi.h \
	ref_cavity.h ref_cell.h ref_cloud.h \
	ref_clump.h ref_collapse.h ref_comprow.h \
	ref_dict.h ref_defs.h ref_edge.h ref_elast.h ref_export.h \
	ref_face.h ref_fixture.h ref_fortran.h \
	ref_gather.h ref_geom.h ref_grid.h \
	ref_histogram.h ref_html.h \
	ref_import.h ref_inflate.h ref_interp.h \
	ref_list.h ref_layer.h \
	ref_malloc.h \
	ref_math.h ref_matrix.h ref_metric.h ref_migrate.h ref_mpi.h \
	ref_node.h ref_part.h ref_phys.h ref_recon.h \
	ref_search.h ref_shard.h ref_smooth.h ref_sort.h ref_split.h \
	ref_subdiv.h ref_swap.h ref_twod.h ref_validation.h 

partioner_ldadd = @zoltan_ldadd@ @parmetis_ldadd@

libref2_a_CFLAGS = @zoltan_include@ @parmetis_include@ @egads_include@
if BUILD_MPI
libref2_a_CFLAGS += -I@MPIINC@
endif # BUILD_MPI

bin_PROGRAMS =

bin_PROGRAMS += ref
ref_SOURCES = ref_subcommand.c
ref_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

bin_PROGRAMS += ref_driver
ref_driver_SOURCES = ref_driver.c
ref_driver_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

bin_PROGRAMS += ref_translate
ref_translate_SOURCES = ref_translate.c
ref_translate_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

bin_PROGRAMS += ref_axisymmetric
ref_axisymmetric_SOURCES = ref_axisymmetric.c
ref_axisymmetric_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

bin_PROGRAMS += ref_inflatable
ref_inflatable_SOURCES = ref_inflatable.c
ref_inflatable_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

bin_PROGRAMS += ref_acceptance
ref_acceptance_SOURCES = ref_acceptance.c
ref_acceptance_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

bin_PROGRAMS += ref_bamg_wrapper
ref_bamg_wrapper_SOURCES = ref_bamg_wrapper.c
ref_bamg_wrapper_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

bin_PROGRAMS += ref_bamg_emulator
ref_bamg_emulator_SOURCES = ref_bamg_emulator.c
ref_bamg_emulator_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

bin_PROGRAMS += ref_fefloa2d
ref_fefloa2d_SOURCES = ref_fefloa2d.c
ref_fefloa2d_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

noinst_PROGRAMS =
TESTS =

TESTS += ref_adapt_test
noinst_PROGRAMS += ref_adapt_test
ref_adapt_test_SOURCES = ref_adapt_test.c
ref_adapt_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_adj_test
noinst_PROGRAMS += ref_adj_test
ref_adj_test_SOURCES = ref_adj_test.c
ref_adj_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_agents_test
noinst_PROGRAMS += ref_agents_test
ref_agents_test_SOURCES = ref_agents_test.c
ref_agents_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_args_test
noinst_PROGRAMS += ref_args_test
ref_args_test_SOURCES = ref_args_test.c
ref_args_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_axi_test
noinst_PROGRAMS += ref_axi_test
ref_axi_test_SOURCES = ref_axi_test.c
ref_axi_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_cavity_test
noinst_PROGRAMS += ref_cavity_test
ref_cavity_test_SOURCES = ref_cavity_test.c
ref_cavity_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_cell_test
noinst_PROGRAMS += ref_cell_test
ref_cell_test_SOURCES = ref_cell_test.c
ref_cell_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_cloud_test
noinst_PROGRAMS += ref_cloud_test
ref_cloud_test_SOURCES = ref_cloud_test.c
ref_cloud_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_clump_test
noinst_PROGRAMS += ref_clump_test
ref_clump_test_SOURCES = ref_clump_test.c
ref_clump_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_collapse_test
noinst_PROGRAMS += ref_collapse_test
ref_collapse_test_SOURCES = ref_collapse_test.c
ref_collapse_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_comprow_test
noinst_PROGRAMS += ref_comprow_test
ref_comprow_test_SOURCES = ref_comprow_test.c
ref_comprow_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_dict_test
noinst_PROGRAMS += ref_dict_test
ref_dict_test_SOURCES = ref_dict_test.c
ref_dict_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_edge_test
noinst_PROGRAMS += ref_edge_test
ref_edge_test_SOURCES = ref_edge_test.c
ref_edge_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_elast_test
noinst_PROGRAMS += ref_elast_test
ref_elast_test_SOURCES = ref_elast_test.c
ref_elast_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_export_test
noinst_PROGRAMS += ref_export_test
ref_export_test_SOURCES = ref_export_test.c
ref_export_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_face_test
noinst_PROGRAMS += ref_face_test
ref_face_test_SOURCES = ref_face_test.c
ref_face_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_fixture_test
noinst_PROGRAMS += ref_fixture_test
ref_fixture_test_SOURCES = ref_fixture_test.c
ref_fixture_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_fortran_test
noinst_PROGRAMS += ref_fortran_test
ref_fortran_test_SOURCES = ref_fortran_test.c
ref_fortran_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_gather_test
noinst_PROGRAMS += ref_gather_test
ref_gather_test_SOURCES = ref_gather_test.c
ref_gather_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_geom_test
noinst_PROGRAMS += ref_geom_test
ref_geom_test_SOURCES = ref_geom_test.c
ref_geom_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_grid_test
noinst_PROGRAMS += ref_grid_test
ref_grid_test_SOURCES = ref_grid_test.c
ref_grid_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_histogram_test
noinst_PROGRAMS += ref_histogram_test
ref_histogram_test_SOURCES = ref_histogram_test.c
ref_histogram_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_html_test
noinst_PROGRAMS += ref_html_test
ref_html_test_SOURCES = ref_html_test.c
ref_html_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_import_test
noinst_PROGRAMS += ref_import_test
ref_import_test_SOURCES = ref_import_test.c
ref_import_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_inflate_test
noinst_PROGRAMS += ref_inflate_test
ref_inflate_test_SOURCES = ref_inflate_test.c
ref_inflate_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_interp_test
noinst_PROGRAMS += ref_interp_test
ref_interp_test_SOURCES = ref_interp_test.c
ref_interp_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_layer_test
noinst_PROGRAMS += ref_layer_test
ref_layer_test_SOURCES = ref_layer_test.c
ref_layer_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_list_test
noinst_PROGRAMS += ref_list_test
ref_list_test_SOURCES = ref_list_test.c
ref_list_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_math_test
noinst_PROGRAMS += ref_math_test
ref_math_test_SOURCES = ref_math_test.c
ref_math_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_matrix_test
noinst_PROGRAMS += ref_matrix_test
ref_matrix_test_SOURCES = ref_matrix_test.c
ref_matrix_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_metric_test
noinst_PROGRAMS += ref_metric_test
ref_metric_test_SOURCES = ref_metric_test.c
ref_metric_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_migrate_test
noinst_PROGRAMS += ref_migrate_test
ref_migrate_test_SOURCES = ref_migrate_test.c
ref_migrate_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_mpi_test
noinst_PROGRAMS += ref_mpi_test
ref_mpi_test_SOURCES = ref_mpi_test.c
ref_mpi_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_node_test
noinst_PROGRAMS += ref_node_test
ref_node_test_SOURCES = ref_node_test.c
ref_node_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_part_test
noinst_PROGRAMS += ref_part_test
ref_part_test_SOURCES = ref_part_test.c
ref_part_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_phys_test
noinst_PROGRAMS += ref_phys_test
ref_phys_test_SOURCES = ref_phys_test.c
ref_phys_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_recon_test
noinst_PROGRAMS += ref_recon_test
ref_recon_test_SOURCES = ref_recon_test.c
ref_recon_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_search_test
noinst_PROGRAMS += ref_search_test
ref_search_test_SOURCES = ref_search_test.c
ref_search_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_shard_test
noinst_PROGRAMS += ref_shard_test
ref_shard_test_SOURCES = ref_shard_test.c
ref_shard_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_smooth_test
noinst_PROGRAMS += ref_smooth_test
ref_smooth_test_SOURCES = ref_smooth_test.c
ref_smooth_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_sort_test
noinst_PROGRAMS += ref_sort_test
ref_sort_test_SOURCES = ref_sort_test.c
ref_sort_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_split_test
noinst_PROGRAMS += ref_split_test
ref_split_test_SOURCES = ref_split_test.c
ref_split_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_subdiv_test
noinst_PROGRAMS += ref_subdiv_test
ref_subdiv_test_SOURCES = ref_subdiv_test.c
ref_subdiv_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ $(partioner_ldadd) @mpi_ldadd@ -lm

TESTS += ref_swap_test
noinst_PROGRAMS += ref_swap_test
ref_swap_test_SOURCES = ref_swap_test.c
ref_swap_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_twod_test
noinst_PROGRAMS += ref_twod_test
ref_twod_test_SOURCES = ref_twod_test.c
ref_twod_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm

TESTS += ref_validation_test
noinst_PROGRAMS += ref_validation_test
ref_validation_test_SOURCES = ref_validation_test.c
ref_validation_test_LDADD = libref2.a @egads_ldadd@ @opencascade_ldadd@ @mpi_ldadd@ -lm
