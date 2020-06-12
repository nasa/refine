function(set_target_rpath NAME)
    if(APPLE)
        set(CMAKE_MACOSX_RPATH TRUE)
    endif()

    set_target_properties(${NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)

    # Don't skip the full RPATH for the build tree
    set_target_properties(${NAME} PROPERTIES SKIP_BUILD_RPATH FALSE)

    # add the automatically determined parts of the RPATH
    set_target_properties(${NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

    foreach(rpath_link_dir ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
        # append RPATH to be used when installing, but only if it's not a system directory
        list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${rpath_link_dir}" isSystemDir)
        if("${isSystemDir}" STREQUAL "-1")
            list(APPEND target_rpath_directories ${rpath_link_dir})
        endif("${isSystemDir}" STREQUAL "-1")
    endforeach()
    set_target_properties(${NAME} PROPERTIES INSTALL_RPATH "${target_rpath_directories}")
endfunction()
