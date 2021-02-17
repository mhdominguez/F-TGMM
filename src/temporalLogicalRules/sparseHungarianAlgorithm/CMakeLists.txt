file(GLOB SRC_HDRS *.h)
file(GLOB SRC_C *.c)
file(GLOB SRC_CPP *.cpp)

file(GLOB_RECURSE ExternalC_SRCS external/*.c)
file(GLOB_RECURSE ExternalCXX_SRCS external/*.cpp)
file(GLOB_RECURSE External_HDRS external/*.h)

set(tgt SparseHungarianAlgorithm)
add_library(${tgt} 
    ${SRC_CPP} ${SRC_C}
    ${SRC_HDRS} ${ExternalC_SRCS} 
    ${ExternalCXX_SRCS} ${ExternalC_HDRS})
target_compile_definitions(${tgt} PRIVATE IL_STD)
target_include_directories(${tgt} PRIVATE external)
set_property(TARGET ${tgt} PROPERTY CXX_STANDARD 11)

IF(WIN32)
    target_compile_options(${tgt} PRIVATE "-W4")
ELSE()
    target_compile_options(${tgt} PRIVATE "-Wall") #?necessary=> -std=gnu++0x")
ENDIF()

set(SPARSEHUNGARIAN_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)


#organize IDE
source_group("Header Files"  FILES ${SRC_HDRS})
source_group("Source Files"  FILES ${SRC_CPP})
source_group("Source Files"  FILES ${SRC_C})


# utilitiy function
function(target_add_sparseHungarianAlgorithm tgt)
    target_include_directories(${tgt} PUBLIC ${SPARSEHUNGARIAN_INCLUDE_DIR})
    target_link_libraries(${tgt} SparseHungarianAlgorithm)
endfunction()
