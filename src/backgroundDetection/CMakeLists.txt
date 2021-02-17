file(GLOB SRC_HDRS *.h)
file(GLOB SRC_CPP *.cpp)

# targets 
set(tgt BackgroundDetectionLib)
add_library(${tgt} ${SRC_CPP} ${SRC_HDRS} )
target_add_mylib(${tgt})
target_add_temporalLogicalRules(${tgt})
target_add_ellipticalHaarFeatures(${tgt})

# utilities
set(BACKGROUNDDETECTION_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

function(target_add_BackgroundDetectionLib tgt)
    target_include_directories(${tgt} PUBLIC ${BACKGROUNDDETECTION_INCLUDE_DIR})
    target_link_libraries(${tgt} BackgroundDetectionLib)
endfunction()
