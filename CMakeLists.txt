
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(refine)

option(refine_ENABLE_TESTS "build refine tests and enable ctest support" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(RPathHelpers)

option(ENABLE_OpenCASCADE "Enable OpenCASCADE support" ON)
if(ENABLE_OpenCASCADE)
    find_package(OpenCASCADE QUIET
        COMPONENTS FoundationClasses ModelingData ModelingAlgorithms ApplicationFramework DataExchange)
endif()

set(THIRD_PARTY_PACKAGES Zoltan ParMETIS EGADS meshlink MPI)
foreach(TPL ${THIRD_PARTY_PACKAGES})
    option(ENABLE_${TPL} "Enable ${TPL} support" ON)
    if(ENABLE_${TPL})
        find_package(${TPL} QUIET)
    endif()
endforeach()

add_subdirectory(src)

