
cmake_minimum_required(VERSION 3.0)

project(refine)

option(refine_ENABLE_TESTS "build refine tests and enable ctest support" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(RPathHelpers)

set(THIRD_PARTY_PACKAGES Zoltan ParMETIS OpenCASCADE EGADS meshlink MPI)
foreach(TPL ${THIRD_PARTY_PACKAGES})
    option(ENABLE_${TPL} "Enable ${TPL} support" ON)
    if(ENABLE_${TPL})
        find_package(${TPL} QUIET)
    endif()
endforeach()

add_subdirectory(src)

