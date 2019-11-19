
set(REFINE_SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

set(REF_TWO_HEADERS
        ${REFINE_SRC_PATH}/ref_adapt.h
        ${REFINE_SRC_PATH}/ref_adj.h
        ${REFINE_SRC_PATH}/ref_agents.h
        ${REFINE_SRC_PATH}/ref_args.h
        ${REFINE_SRC_PATH}/ref_axi.h
        ${REFINE_SRC_PATH}/ref_cavity.h
        ${REFINE_SRC_PATH}/ref_cell.h
        ${REFINE_SRC_PATH}/ref_cloud.h
        ${REFINE_SRC_PATH}/ref_clump.h
        ${REFINE_SRC_PATH}/ref_collapse.h
        ${REFINE_SRC_PATH}/ref_comprow.h
        ${REFINE_SRC_PATH}/ref_defs.h
        ${REFINE_SRC_PATH}/ref_dict.h
        ${REFINE_SRC_PATH}/ref_edge.h
        ${REFINE_SRC_PATH}/ref_elast.h
        ${REFINE_SRC_PATH}/ref_endian.h
        ${REFINE_SRC_PATH}/ref_export.h
        ${REFINE_SRC_PATH}/ref_face.h
        ${REFINE_SRC_PATH}/ref_fixture.h
        ${REFINE_SRC_PATH}/ref_fortran.h
        ${REFINE_SRC_PATH}/ref_gather.h
        ${REFINE_SRC_PATH}/ref_geom.h
        ${REFINE_SRC_PATH}/ref_grid.h
        ${REFINE_SRC_PATH}/ref_histogram.h
        ${REFINE_SRC_PATH}/ref_html.h
        ${REFINE_SRC_PATH}/ref_import.h
        ${REFINE_SRC_PATH}/ref_inflate.h
        ${REFINE_SRC_PATH}/ref_interp.h
        ${REFINE_SRC_PATH}/ref_layer.h
        ${REFINE_SRC_PATH}/ref_list.h
        ${REFINE_SRC_PATH}/ref_malloc.h
        ${REFINE_SRC_PATH}/ref_math.h
        ${REFINE_SRC_PATH}/ref_matrix.h
        ${REFINE_SRC_PATH}/ref_metric.h
        ${REFINE_SRC_PATH}/ref_migrate.h
        ${REFINE_SRC_PATH}/ref_mpi.h
        ${REFINE_SRC_PATH}/ref_node.h
        ${REFINE_SRC_PATH}/ref_part.h
        ${REFINE_SRC_PATH}/ref_phys.h
        ${REFINE_SRC_PATH}/ref_recon.h
        ${REFINE_SRC_PATH}/ref_search.h
        ${REFINE_SRC_PATH}/ref_shard.h
        ${REFINE_SRC_PATH}/ref_smooth.h
        ${REFINE_SRC_PATH}/ref_sort.h
        ${REFINE_SRC_PATH}/ref_split.h
        ${REFINE_SRC_PATH}/ref_subdiv.h
        ${REFINE_SRC_PATH}/ref_swap.h
        ${REFINE_SRC_PATH}/ref_validation.h
        )

set(REF_TWO_SRC
  "${REFINE_SRC_PATH}/ref_adapt.c"
  "${REFINE_SRC_PATH}/ref_agents.c"
  "${REFINE_SRC_PATH}/ref_adj.c"
  "${REFINE_SRC_PATH}/ref_args.c"
  "${REFINE_SRC_PATH}/ref_axi.c"
  "${REFINE_SRC_PATH}/ref_cavity.c"
  "${REFINE_SRC_PATH}/ref_cell.c"
  "${REFINE_SRC_PATH}/ref_cloud.c"
  "${REFINE_SRC_PATH}/ref_clump.c"
  "${REFINE_SRC_PATH}/ref_collapse.c"
  "${REFINE_SRC_PATH}/ref_comprow.c"
  "${REFINE_SRC_PATH}/ref_dict.c"
  "${REFINE_SRC_PATH}/ref_edge.c"
  "${REFINE_SRC_PATH}/ref_elast.h"
  "${REFINE_SRC_PATH}/ref_endian.h"
  "${REFINE_SRC_PATH}/ref_export.c"
  "${REFINE_SRC_PATH}/ref_face.c"
  "${REFINE_SRC_PATH}/ref_fixture.c"
  "${REFINE_SRC_PATH}/ref_fortran.c"
  "${REFINE_SRC_PATH}/ref_gather.c"
  "${REFINE_SRC_PATH}/ref_geom.c"
  "${REFINE_SRC_PATH}/ref_grid.c"
  "${REFINE_SRC_PATH}/ref_histogram.c"
  "${REFINE_SRC_PATH}/ref_html.c"
  "${REFINE_SRC_PATH}/ref_import.c"
  "${REFINE_SRC_PATH}/ref_inflate.c"
  "${REFINE_SRC_PATH}/ref_interp.c"
  "${REFINE_SRC_PATH}/ref_list.c"
  "${REFINE_SRC_PATH}/ref_math.c"
  "${REFINE_SRC_PATH}/ref_matrix.c"
  "${REFINE_SRC_PATH}/ref_metric.c"
  "${REFINE_SRC_PATH}/ref_migrate.c"
  "${REFINE_SRC_PATH}/ref_mpi.c"
  "${REFINE_SRC_PATH}/ref_node.c"
  "${REFINE_SRC_PATH}/ref_part.c"
  "${REFINE_SRC_PATH}/ref_phys.c"
  "${REFINE_SRC_PATH}/ref_recon.c"
  "${REFINE_SRC_PATH}/ref_search.c"
  "${REFINE_SRC_PATH}/ref_shard.c"
  "${REFINE_SRC_PATH}/ref_smooth.c"
  "${REFINE_SRC_PATH}/ref_sort.c"
  "${REFINE_SRC_PATH}/ref_split.c"
  "${REFINE_SRC_PATH}/ref_subdiv.c"
  "${REFINE_SRC_PATH}/ref_swap.c"
  "${REFINE_SRC_PATH}/ref_validation.c"
)

add_executable(ref_phys_test ${REF_TWO_SRC} ${REF_TWO_HEADERS} ${REFINE_SRC_PATH}/ref_phys_test.c)
add_executable(ref_migrate_test ${REF_TWO_SRC} ${REF_TWO_HEADERS} ${REFINE_SRC_PATH}/ref_migrate_test.c)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ref_phys_test PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
  target_compile_options(ref_migrate_test PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
endif()
