set(REF_TWO_HEADERS
        ref_adapt.h
        ref_adj.h
        ref_agents.h
        ref_args.h
        ref_axi.h
        ref_cavity.h
        ref_cell.h
        ref_cloud.h
        ref_clump.h
        ref_collapse.h
        ref_comprow.h
        ref_defs.h
        ref_dict.h
        ref_dist.h
        ref_edge.h
        ref_egads.h
        ref_elast.h
        ref_endian.h
        ref_export.h
        ref_face.h
        ref_fixture.h
        ref_fortran.h
        ref_gather.h
        ref_geom.h
        ref_grid.h
        ref_histogram.h
        ref_html.h
        ref_import.h
        ref_inflate.h
        ref_interp.h
        ref_layer.h
        ref_list.h
        ref_malloc.h
        ref_math.h
        ref_matrix.h
        ref_meshlink.h
        ref_metric.h
        ref_migrate.h
        ref_mpi.h
        ref_node.h
        ref_part.h
        ref_phys.h
        ref_recon.h
        ref_search.h
        ref_shard.h
        ref_smooth.h
        ref_sort.h
        ref_split.h
        ref_subdiv.h
        ref_swap.h
        ref_validation.h
        ref_elast.h
        ref_endian.h
        )

set(REF_CORE_SRC
        ref_adapt.c
        ref_agents.c
        ref_adj.c
        ref_args.c
        ref_axi.c
        ref_cavity.c
        ref_cell.c
        ref_cloud.c
        ref_clump.c
        ref_collapse.c
        ref_comprow.c
        ref_dict.c
        ref_dist.c
        ref_edge.c
        ref_export.c
        ref_face.c
        ref_fixture.c
        ref_fortran.c
        ref_gather.c
        ref_geom.c
        ref_grid.c
        ref_histogram.c
        ref_html.c
        ref_import.c
        ref_inflate.c
        ref_interp.c
        ref_list.c
        ref_math.c
        ref_matrix.c
        ref_meshlink.c
        ref_metric.c
        ref_node.c
        ref_part.c
        ref_phys.c
        ref_recon.c
        ref_search.c
        ref_shard.c
        ref_smooth.c
        ref_sort.c
        ref_split.c
        ref_subdiv.c
        ref_swap.c
        ref_validation.c
        )

function(install_target TARGET_NAME)
    install(TARGETS ${TARGET_NAME}
            DESTINATION bin)
    set_target_rpath(${TARGET_NAME} ../lib)
endfunction(install_target)

function(create_program TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${TARGET_NAME} PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
    endif ()
    install_target(${TARGET_NAME})
endfunction(create_program)

if (EGADS_FOUND)
    message(STATUS "EGADS Found: ${EGADS_LIBRARIES}")
    message(STATUS "EGADSLITE Found: ${EGADSLITE_LIBRARIES}")
    list(APPEND THIRD_PARTY_INCLUDES ${EGADS_INCLUDE_DIR})
    list(APPEND EXTRA_DEFINITIONS HAVE_EGADS)
endif ()

if (PARMETIS_FOUND)
    message(STATUS "ParMETIS Found: ${PARMETIS_LIBRARIES}")
    list(APPEND THIRD_PARTY_LIBRARIES ${PARMETIS_LIBRARIES})
    list(APPEND THIRD_PARTY_INCLUDES ${PARMETIS_INCLUDE_DIRS})
    list(APPEND EXTRA_DEFINITIONS HAVE_PARMETIS)
endif ()

if (ZOLTAN_FOUND)
    message(STATUS "Zoltan Found: ${ZOLTAN_LIBRARIES}")
    list(APPEND THIRD_PARTY_LIBRARIES ${ZOLTAN_LIBRARIES})
    list(APPEND THIRD_PARTY_INCLUDES ${ZOLTAN_INCLUDE_DIRS})
    list(APPEND EXTRA_DEFINITIONS HAVE_ZOLTAN)
endif ()

find_library(MATH_LIBRARY m REQUIRED)
if (MATH_LIBRARY)
    list(APPEND THIRD_PARTY_LIBRARIES ${MATH_LIBRARY})
endif ()

set(REF_MPI_SRC
        ref_mpi.c
        ref_migrate.c
        )

add_library(refine_with_egadslite STATIC ref_egads.c)
target_compile_definitions(refine_with_egadslite PRIVATE HAVE_EGADS_LITE)

add_library(refine_with_egads STATIC ref_egads.c)
if (EGADS_FOUND)
    target_link_libraries(refine_with_egads PRIVATE ${EGADS_LIBRARIES})
    target_include_directories(refine_with_egads PRIVATE ${EGADS_INCLUDE_DIRS})
    target_compile_definitions(refine_with_egads PRIVATE HAVE_EGADS)
    if (OpenCASCADE_FOUND)
        message(STATUS "OpenCASCADE Found: ${OpenCASCADE_LIBRARIES}")
        target_link_libraries(refine_with_egads PRIVATE ${OpenCASCADE_LIBRARIES})
    endif ()

    target_compile_definitions(refine_with_egadslite PRIVATE HAVE_EGADS)
    target_include_directories(refine_with_egadslite PRIVATE ${EGADS_INCLUDE_DIRS})
    target_link_libraries(refine_with_egadslite PRIVATE ${EGADSLITE_LIBRARIES})
endif ()

add_library(refine_without_mpi STATIC ${REF_MPI_SRC})
if (MPI_FOUND)
    add_library(refine_with_mpi STATIC ${REF_MPI_SRC})
    target_include_directories(refine_with_mpi PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries(refine_with_mpi PRIVATE ${MPI_LIBRARIES})
    target_compile_definitions(refine_with_mpi PRIVATE HAVE_MPI)
endif ()

add_library(refine_core STATIC ${REF_CORE_SRC} ${REF_TWO_HEADERS})
target_link_libraries(refine_core PRIVATE ${THIRD_PARTY_LIBRARIES})
target_include_directories(refine_core PRIVATE ${THIRD_PARTY_INCLUDES})
target_compile_definitions(refine_core PRIVATE ${EXTRA_DEFINITIONS})

create_program(ref ref_subcommand.c)
target_link_libraries(ref PRIVATE refine_core refine_without_mpi refine_with_egads)

if (MPI_FOUND)
    message("-- MPI Found: ${MPI_LIBRARIES}")
    message("-- building refmpifull with MPI support")
    create_program(refmpifull ref_subcommand.c)
    target_link_libraries(refmpifull PRIVATE refine_core refine_with_mpi refine_with_egads)

    message("-- building refmpi with MPI support")
    create_program(refmpi ref_subcommand.c)
    target_link_libraries(refmpi PRIVATE refine_core refine_with_mpi refine_with_egadslite)
endif ()
