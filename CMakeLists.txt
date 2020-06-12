cmake_minimum_required(VERSION 3.0)

project(refine)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(RPathHelpers)

find_package(MPI QUIET)
find_package(OpenCASCADE QUIET)
find_package(PARMETIS QUIET)
find_package(Zoltan QUIET)
find_package(EGADS QUIET)

add_subdirectory(src)
