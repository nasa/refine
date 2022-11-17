cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(refine)

find_program(GIT_EXEC git)
if (GIT_EXEC AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    execute_process(COMMAND sh -c "${GIT_EXEC} describe --tag | tr -d '\n'"
            OUTPUT_VARIABLE REFINE_VERSION)
    message(STATUS "refine version: ${REFINE_VERSION}")
else()
    set(REFINE_VERSION 3.0.0)
    message(STATUS "not git repo, guess refine version: ${REFINE_VERSION}")
endif()

option(refine_ENABLE_TESTS "build refine tests and enable ctest support" ON)
option(refine_ENABLE_CUDA "build refine with CUDA support (requires nvcc compiler)" OFF)

if (refine_ENABLE_CUDA)
    cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    option(CMAKE_CUDA_SEPARABLE_COMPILATION "Use separable compilation" ON)
    # option(CMAKE_CUDA_RESOLVE_DEVICE_SYMBOLS ON CACHE "" BOOL)
    option(CMAKE_CUDA_ARCHITECTURES "Default to NVIDIA VOLTA GPU support" 70)
endif()

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
