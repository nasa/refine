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
	set_target_rpath(${TARGET_NAME} ../lib)
endfunction(install_target)

function (create_program TARGET_NAME REFINE_LIBRARY)
	add_executable(${TARGET_NAME} ${ARGN})
	target_link_libraries(${TARGET_NAME} PRIVATE ${REFINE_LIBRARY})
	if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
		target_compile_options(${TARGET_NAME} PRIVATE -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized)
	endif()
	install_target(${TARGET_NAME})
endfunction(create_program)

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
target_link_libraries(refine PRIVATE ${THIRD_PARTY_LIBRARIES})
target_include_directories(refine PRIVATE ${THIRD_PARTY_INCLUDES})
target_compile_definitions(refine PRIVATE ${EXTRA_DEFINITIONS})

set(PROGRAMS_THAT_DO_NOT_REQUIRE_MPI
		ref_driver
		ref_translate
		ref_axisymmetric
		ref_inflatable
		ref_acceptance
		)

foreach(PROGRAM ${PROGRAMS_THAT_DO_NOT_REQUIRE_MPI})
	create_program(${PROGRAM} refine ${PROGRAM}.c)
endforeach()
create_program(ref refine ref_subcommand.c)

if (MPI_FOUND)
	message("-- MPI Found: ${MPI_LIBRARIES}")

	add_library(refinempi STATIC ${REF_TWO_SRC} ${REF_TWO_HEADERS})
	target_link_libraries(refinempi PRIVATE ${THIRD_PARTY_LIBRARIES} ${MPI_LIBRARIES})
	target_include_directories(refinempi PRIVATE ${THIRD_PARTY_INCLUDES} ${MPI_INCLUDE_PATH})
	target_compile_definitions(refinempi PRIVATE ${EXTRA_DEFINITIONS} HAVE_MPI)

	message("-- building refmpifull with MPI support")
	create_program(refmpifull refinempi ref_subcommand.c)
endif()
