file(GLOB PROJ_HDRS *.h)
file(GLOB PROJ_SRCS *.cpp)

add_library(MedianFilter2D_lib ${PROJ_SRCS} ${PROJ_HDRS})

# utlitiles
set(MEDIANFILTER2D_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

function(target_add_MedianFilter2D tgt)
    target_include_directories(${tgt} PUBLIC ${MEDIANFILTER2D_INCLUDE_DIR})
    target_link_libraries(${tgt} MedianFilter2D_lib)
endfunction()