find_path(ZOLTAN_INCLUDE_DIR zoltan.h)
find_library(ZOLTAN_LIBRARY zoltan)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
        ZOLTAN
        DEFAULT_MSG
        ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR
)

mark_as_advanced(ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY)

if(ZOLTAN_FOUND AND NOT TARGET ZOLTAN::ZOLTAN)
    message(STATUS "ZOLTAN Found: ${ZOLTAN_LIBRARY}")
    add_library(ZOLTAN::ZOLTAN UNKNOWN IMPORTED)
    set_target_properties(ZOLTAN::ZOLTAN PROPERTIES
            IMPORTED_LOCATION ${ZOLTAN_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${ZOLTAN_INCLUDE_DIR}
            )
endif()
