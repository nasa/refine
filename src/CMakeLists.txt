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

set(REF_TWO_SRC
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
    ref_egads.c
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
    ref_migrate.c
    ref_mpi.c
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
    set_target_rpath(${TARGET_NAME} ${CMAKE_INSTALL_PREFIX})
endfunction()

if (EGADS_FOUND)
	message(STATUS "EGADS Found: ${EGADS_LIBRARIES}")
	list(APPEND THIRD_PARTY_LIBRARIES ${EGADS_LIBRARIES})
	list(APPEND THIRD_PARTY_INCLUDES ${EGADS_INCLUDE_DIR})
	list(APPEND EXTRA_DEFINITIONS HAVE_EGADS)
endif()

if (PARMETIS_FOUND)
	message(STATUS "ParMETIS Found: ${PARMETIS_LIBRARIES}")
	list(APPEND THIRD_PARTY_LIBRARIES ${PARMETIS_LIBRARIES})
	list(APPEND THIRD_PARTY_INCLUDES ${PARMETIS_INCLUDE_DIRS})
	list(APPEND EXTRA_DEFINITIONS HAVE_PARMETIS)
endif()

if (ZOLTAN_FOUND)
	message(STATUS "Zoltan Found: ${ZOLTAN_LIBRARIES}")
	list(APPEND THIRD_PARTY_LIBRARIES ${ZOLTAN_LIBRARIES})
	list(APPEND THIRD_PARTY_INCLUDES ${ZOLTAN_INCLUDE_DIRS})
	list(APPEND EXTRA_DEFINITIONS HAVE_ZOLTAN)
endif()

if (OpenCASCADE_FOUND)
	message(STATUS "OpenCASCADE Found: ${OpenCASCADE_LIBRARIES}")
	list(APPEND THIRD_PARTY_LIBRARIES ${OpenCASCADE_LIBRARIES})
endif()

find_library(MATH_LIBRARY m REQUIRED)
if(MATH_LIBRARY)
    list(APPEND THIRD_PARTY_LIBRARIES ${MATH_LIBRARY})
endif()

add_library(refine STATIC ${REF_TWO_SRC} ${REF_TWO_HEADERS})
target_link_libraries(refine PUBLIC ${THIRD_PARTY_LIBRARIES})
target_include_directories(refine PUBLIC ${THIRD_PARTY_INCLUDES})
target_compile_definitions(refine PUBLIC ${EXTRA_DEFINITIONS})

add_executable(ref ref_subcommand.c)
target_link_libraries(ref PRIVATE refine)

add_executable(ref_phys_test ref_phys_test.c)
target_link_libraries(ref_phys_test PRIVATE refine)

add_executable(ref_migrate_test ref_migrate_test.c)
target_link_libraries(ref_migrate_test PRIVATE refine)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ref_phys_test PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
  target_compile_options(ref_migrate_test PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
endif()

install_target(ref)
install_target(ref_phys_test)
install_target(ref_migrate_test)

if (MPI_FOUND)
	message("-- MPI Found: ${MPI_LIBRARIES}")
	message("-- building refmpifull with MPI support")
	add_executable(refmpifull ${REF_TWO_SRC} ${REF_TWO_HEADERS} ref_subcommand.c)
	target_link_libraries(refmpifull PRIVATE ${THIRD_PARTY_LIBRARIES} ${MPI_LIBRARIES})
	target_compile_definitions(refmpifull PRIVATE ${EXTRA_DEFINITIONS} HAVE_MPI)
	target_include_directories(refmpifull PRIVATE ${THIRD_PARTY_INCLUDES} ${MPI_INCLUDE_PATH})
    install_target(refmpifull)
endif()
