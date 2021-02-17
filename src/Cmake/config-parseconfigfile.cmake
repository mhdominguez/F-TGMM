add_library(parse_config_file 
    Utils/parseConfigFile.cpp
    Utils/parseConfigFile.h)

set(PARSECONFIGFILE_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/../Utils")

function(target_add_parseconfigfile tgt)
    target_include_directories(${tgt} PUBLIC ${PARSECONFIGFILE_INCLUDE_DIR})
    target_link_libraries(${tgt} parse_config_file)
endfunction()