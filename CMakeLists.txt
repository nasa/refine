cmake_minimum_required(VERSION 3.0)

project(refine)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(OpenCASCADE REQUIRED)
find_package(EGADS REQUIRED)

add_subdirectory(src)
