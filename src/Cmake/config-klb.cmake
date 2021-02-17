add_subdirectory(keller-lab-block-filetype/src)
set(KLB_INCLUDE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../keller-lab-block-filetype/src)

function(target_add_klb tgt)
    target_include_directories(${tgt} PUBLIC ${KLB_INCLUDE_DIRECTORY})
    target_link_libraries(${tgt} klb_static)
endfunction()
