
cmake_minimum_required(VERSION 3.0)

project(refine)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(RPathHelpers)
include(CTest)
enable_testing()

set(THIRD_PARTY_PACKAGES Zoltan ParMETIS OpenCASCADE EGADS MPI)
foreach(TPL ${THIRD_PARTY_PACKAGES})
    option(ENABLE_${TPL} "Enable ${TPL} support" ON)
    if (ENABLE_${TPL})
        find_package(${TPL} QUIET)
    endif()
endforeach()

add_subdirectory(src)

