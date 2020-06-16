# - Try to find parmetis
# Once done this will define
#  PARMETIS_FOUND - System has PARMETIS
#  PARMETIS_INCLUDE_DIRS - The PARMETIS include directories
#  PARMETIS_LIBRARIES - The libraries needed to use PARMETIS
#  PARMETIS_DEFINITIONS - Compiler switches required for using PARMETIS

find_package(METIS)
if(METIS_FOUND)
    find_path(PARMETIS_INCLUDE_DIR parmetis.h)
    find_library(PARMETIS_LIBRARY parmetis)

    set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARIES})
    set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIRS})

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(
            PARMETIS
            DEFAULT_MSG
            PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR
    )

    mark_as_advanced(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY )
else()
    set(PARMETIS_FOUND FALSE)
endif()


