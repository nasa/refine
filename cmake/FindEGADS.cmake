# - Try to find egads
# Once done this will define
#  EGADS_FOUND - System has EGADS
#  EGADS_INCLUDE_DIRS - The EGADS include directories
#  EGADS_LIBRARIES - The libraries needed to use EGADS
#  EGADS_DEFINITIONS - Compiler switches required for using EGADS

find_path(EGADS_INCLUDE_DIR egads.h)

find_library(EGADS_LIBRARY egads)
find_library(EGADSLITE_LIBRARY egadslite)

set(EGADS_LIBRARIES ${EGADS_LIBRARY} )
set(EGADSLITE_LIBRARIES ${EGADSLITE_LIBRARY} )
set(EGADS_INCLUDE_DIRS ${EGADS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EGADS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
        EGADS
        DEFAULT_MSG
        EGADS_LIBRARY EGADSLITE_LIBRARY EGADS_INCLUDE_DIR
)

mark_as_advanced(EGADS_INCLUDE_DIR EGADS_LIBRARY EGADSLITE_LIBRARY)

