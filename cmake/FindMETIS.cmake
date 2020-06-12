# - Try to find metis
# Once done this will define
#  METIS_FOUND - System has METIS
#  METIS_INCLUDE_DIRS - The METIS include directories
#  METIS_LIBRARIES - The libraries needed to use METIS
#  METIS_DEFINITIONS - Compiler switches required for using METIS

set(_PREFIX "${METIS_PREFIX_DEFAULT}" CACHE STRING "Zoltan install directory")
if(METIS_PREFIX)
    message(STATUS "METIS_PREFIX ${METIS_PREFIX}")
endif()

find_path(METIS_INCLUDE_DIR metis.h PATHS "${METIS_PREFIX}/include")

find_library(METIS_LIBRARY metis PATHS "${METIS_PREFIX}/lib")

set(METIS_LIBRARIES ${METIS_LIBRARY} )
set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set METIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
        METIS
        DEFAULT_MSG
        METIS_LIBRARY METIS_INCLUDE_DIR
)

mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY )

