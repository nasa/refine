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
            ParMETIS
            DEFAULT_MSG
            PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR
    )

    mark_as_advanced(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY)
    if(ParMETIS_FOUND AND NOT TARGET ParMETIS::ParMETIS)
        message(STATUS "PARMETIS Found: ${PARMETIS_LIBRARY}")
        add_library(ParMETIS::ParMETIS UNKNOWN IMPORTED)
        set_target_properties(ParMETIS::ParMETIS PROPERTIES
                IMPORTED_LOCATION ${PARMETIS_LIBRARY}
                INTERFACE_INCLUDE_DIRECTORIES ${PARMETIS_INCLUDE_DIR}
                )
    endif()
else()
    set(ParMETIS_FOUND FALSE)
endif()


