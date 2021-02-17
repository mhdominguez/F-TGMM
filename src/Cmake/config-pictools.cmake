#for JP2000 format using PicTools

# Note: This is completely untested.
# 
#   I wrote this during a re-write of the old cmake build after pictools support
#   had been mostly removed.
# 
#   The basics are here, and I think it should be straightforward to get it right
#   if pictools needs to be added back.
# 
#   Nathan Clack, Jan 2016

set(HAVE_PICTOOLS false)
IF( IS_DIRECTORY PICTools )
    add_subdirectory(PicTools)
    set(HAVE_PICTOOLS true)
    MESSAGE("Compiling TGMM code with PicTools for JP2 support. Include dirs ${PICTOOLS_INCLUDE_DIR}")

    #copy DLLs for PICTOOLS
    file(GLOB PicTools_ConfigFiles "PICTools/lib/picx*.dll")
    foreach(ConfigFile ${ConfigFiles})
        install(FILES ${ConfigFile} DESTINATION bin)
    endforeach()

ELSE()
    MESSAGE("Compiling TGMM code WITHOUT PicTools for JP2 support. Folder PICTools not found")
ENDIF()


function(target_add_pictools tgt)
    if(HAVE_PICTOOLS)
        target_link_libraries(${tgt} PicToolsJ2K)
        target_compile_definitions(${tgt} -DPICTOOLS_JP2K)
        target_include_directories(${tgt} ${PICTOOLS_INCLUDE_DIR})

        #copy DLLs for PICTOOLS to build dir
        file(GLOB ConfigFiles "PICTools/lib/picx*.dll")
        foreach(ConfigFile ${ConfigFiles})
            add_custom_command(TARGET ${tgt} POST_BUILD COMMAND 
                ${CMAKE_COMMAND} -E copy ${ConfigFile} $<TARGET_FILE_DIR:${tgt}>)
        endforeach()

    endif()
endfunction()


# OLD STUFF

    #include_directories(PICTools)
    #add_subdirectory("${PROJECT_SOURCE_DIR}/../PICTools" "${CMAKE_CURRENT_BINARY_DIR}/PICTools")
    #INCLUDE_DIRECTORIES( ${PICTOOLS_INCLUDE_DIR} )
    #LINK_DIRECTORIES(${PICTOOLS_LIBRARY_DIR})