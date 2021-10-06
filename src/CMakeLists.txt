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
        ref_facelift.h
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
        ref_iso.h
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
        ref_elast.c
        ref_export.c
        ref_face.c
        ref_facelift.c
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
        ref_iso.c
        ref_layer.c
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

function(create_program TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${TARGET_NAME} PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
    endif ()
    set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE CXX)
    set_target_rpath(${TARGET_NAME})
    install(TARGETS ${TARGET_NAME} DESTINATION bin)
endfunction()

function(create_library TARGET_NAME)
    add_library(${TARGET_NAME} ${ARGN})
    add_library(refine::${TARGET_NAME} ALIAS ${TARGET_NAME})
    target_include_directories(${TARGET_NAME} PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:include/refine>
            )
    install(TARGETS ${TARGET_NAME}
            EXPORT refine
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            )
endfunction()

if (ParMETIS_FOUND)
    list(APPEND THIRD_PARTY_LIBRARIES ParMETIS::ParMETIS)
    list(APPEND EXTRA_DEFINITIONS HAVE_PARMETIS)
endif ()

if (Zoltan_FOUND)
    list(APPEND THIRD_PARTY_LIBRARIES Zoltan::Zoltan)
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

create_library(refine_with_egadslite STATIC ref_egads.c)
# HAVE_EGADS_LITE nested in HAVE_EGADS block, safe to set
target_compile_definitions(refine_with_egadslite PRIVATE HAVE_EGADS_LITE)

create_library(refine_with_egads STATIC ref_egads.c)
if (EGADS_FOUND)
    target_link_libraries(refine_with_egads PRIVATE EGADS::EGADS)
    target_compile_definitions(refine_with_egads PRIVATE HAVE_EGADS)
    if (OpenCASCADE_FOUND)
        message(STATUS "OpenCASCADE Found: ${OpenCASCADE_LIBRARIES}")
        list(FIND OpenCASCADE_LIBRARIES TKernel TKernel_FOUND)
        if (${TKernel_FOUND} EQUAL -1)
            message(WARNING "TKernel target not found in OpenCASCADE libraries: ${OpenCASCADE_LIBRARIES}\nSkipping linking refine with OpenCASCADE")
        else()
            target_link_libraries(refine_with_egads PRIVATE TKernel)
        endif()
    endif ()

    target_compile_definitions(refine_with_egadslite PRIVATE HAVE_EGADS)
    target_link_libraries(refine_with_egadslite PRIVATE EGADS::EGADSLITE)
endif ()

create_library(refine_without_mpi STATIC ${REF_MPI_SRC})
if (MPI_FOUND)
    create_library(refine_with_mpi STATIC ${REF_MPI_SRC})
    target_include_directories(refine_with_mpi PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries(refine_with_mpi PRIVATE ${MPI_C_LIBRARIES})
    target_link_libraries(refine_with_mpi PRIVATE ${THIRD_PARTY_LIBRARIES})
    target_compile_definitions(refine_with_mpi PRIVATE HAVE_MPI ${EXTRA_DEFINITIONS})
endif ()

create_library(refine_core STATIC ${REF_CORE_SRC} ${REF_TWO_HEADERS})
target_link_libraries(refine_core PRIVATE ${THIRD_PARTY_LIBRARIES})
target_include_directories(refine_core PRIVATE ${THIRD_PARTY_INCLUDES})
target_compile_definitions(refine_core PRIVATE ${EXTRA_DEFINITIONS})

if (meshlink_FOUND)
    message(STATUS "meshlink Found")
    target_link_libraries(refine_core PRIVATE
            meshlink::meshlink
            meshlink::mlparser_xerces
            meshlink::mlkernel_geode
            )
    target_compile_definitions(refine_core PRIVATE HAVE_MESHLINK)
else ()
    message(STATUS "meshlink NOT Found")
endif ()

create_program(ref ref_subcommand.c)
target_link_libraries(ref PRIVATE refine_core refine_without_mpi refine_with_egads)

install(EXPORT refine
        NAMESPACE refine::
        DESTINATION cmake)
install(FILES ${REF_TWO_HEADERS}
        DESTINATION include/refine)
install(FILES refineConfig.cmake
        DESTINATION cmake)

if (MPI_FOUND)
    message("-- MPI Found: ${MPI_C_LIBRARIES}")
    message("-- building refmpifull with MPI support")
    create_program(refmpifull ref_subcommand.c)
    target_link_libraries(refmpifull PRIVATE refine_core refine_with_mpi refine_with_egads)

    message("-- building refmpi with MPI support")
    create_program(refmpi ref_subcommand.c)
    target_link_libraries(refmpi PRIVATE refine_core refine_with_mpi refine_with_egadslite)
endif ()

set(REF_TESTS
        ref_adapt_test.c
        ref_adj_test.c
        ref_agents_test.c
        ref_args_test.c
        ref_axi_test.c
        ref_cavity_test.c
        ref_cell_test.c
        ref_cloud_test.c
        ref_clump_test.c
        ref_collapse_test.c
        ref_comprow_test.c
        ref_dict_test.c
        ref_dist_test.c
        ref_edge_test.c
        ref_egads_test.c
        ref_elast_test.c
        ref_export_test.c
        ref_face_test.c
        ref_facelift_test.c
        ref_fixture_test.c
        ref_fortran_test.c
        ref_gather_test.c
        ref_geom_test.c
        ref_grid_test.c
        ref_histogram_test.c
        ref_html_test.c
        ref_import_test.c
        ref_inflate_test.c
        ref_interp_test.c
        ref_iso_test.c
        ref_layer_test.c
        ref_list_test.c
        ref_math_test.c
        ref_matrix_test.c
        ref_meshlink_test.c
        ref_metric_test.c
        ref_migrate_test.c
        ref_mpi_test.c
        ref_node_test.c
        ref_part_test.c
        ref_phys_test.c
        ref_recon_test.c
        ref_search_test.c
        ref_shard_test.c
        ref_smooth_test.c
        ref_sort_test.c
        ref_split_test.c
        ref_subdiv_test.c
        ref_swap_test.c
        ref_validation_test.c
        )

if (refine_ENABLE_TESTS)
    include(CTest)
    enable_testing()
    foreach (REF_TEST_SRC ${REF_TESTS})
        string(REPLACE ".c" "" REF_TEST_NAME ${REF_TEST_SRC})
        add_executable(${REF_TEST_NAME} ${REF_TEST_SRC})
        target_link_libraries(${REF_TEST_NAME} PRIVATE refine_core refine_without_mpi refine_with_egads)
        add_test(NAME ${REF_TEST_NAME} COMMAND ${REF_TEST_NAME})
    endforeach ()
endif ()
