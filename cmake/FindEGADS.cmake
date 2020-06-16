# - Try to find egads
# Once done this will define
#  EGADS_FOUND - System has EGADS
#  EGADS_INCLUDE_DIRS - The EGADS include directories
#  EGADS_LIBRARIES - The libraries needed to use EGADS
#  EGADS_DEFINITIONS - Compiler switches required for using EGADS

find_path(EGADS_INCLUDE_DIR egads.h)

find_library(EGADS_LIBRARY egads)
find_library(EGADSLITE_LIBRARY egadslite)

set(EGADS_INCLUDE_DIRS ${EGADS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EGADS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
        EGADS
        DEFAULT_MSG
        EGADS_LIBRARY EGADSLITE_LIBRARY EGADS_INCLUDE_DIR
)

if(EGADS_FOUND AND NOT TARGET EGADS::EGADS)
    message(STATUS "EGADS Found: ${EGADS_LIBRARY}")
    add_library(EGADS::EGADS UNKNOWN IMPORTED)
    set_target_properties(EGADS::EGADS PROPERTIES
            IMPORTED_LOCATION ${EGADS_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${EGADS_INCLUDE_DIR}
            )
    if(NOT TARGET EGADS::EGADSLITE)
        message(STATUS "EGADSLITE Found: ${EGADSLITE_LIBRARY}")
        add_library(EGADS::EGADSLITE UNKNOWN IMPORTED)
        set_target_properties(EGADS::EGADSLITE PROPERTIES
                IMPORTED_LOCATION ${EGADSLITE_LIBRARY}
                INTERFACE_INCLUDE_DIRECTORIES ${EGADS_INCLUDE_DIR}
                )
    endif()
endif()

mark_as_advanced(EGADS_INCLUDE_DIR EGADS_LIBRARY EGADSLITE_LIBRARY)

