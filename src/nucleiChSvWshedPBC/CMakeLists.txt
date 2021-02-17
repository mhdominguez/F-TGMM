add_subdirectory(MedianFilter2D)
add_subdirectory(CUDAmedianFilter2D)

function(target_add_MedianFilter2D_CUDAlib tgt)
    if(HAVE_CUDA)
        target_link_libraries(${tgt} MedianFilter2D_CUDAlib)
        target_include_directories(${tgt} PUBLIC ${MEDIANFILTER2DCUDA_INCLUDE_DIR})
    else()
        target_add_MedianFilter2D(${tgt})
    endif()
endfunction()

# 
# submodules: add source and header files from subfolders
# 
file(GLOB SRC_HDRS *.h)
file(GLOB SRC_C *.c)

file(GLOB MEDIAN_SRC_HDRS external/medianFilter2D/*.h)
file(GLOB MEDIAN_SRC_C    external/medianFilter2D/*.c)

SET(SRC_CPP 
    watershedPersistanceAgglomeration.cpp
    hierarchicalSegmentation.cpp
    autoThreshold.cpp)

# targets

set(tgt WatershedWithPersistanceAgglomeration_Test)
add_executable(${tgt} 
    mainTest.cpp
    ${SRC_CPP} ${SRC_HDRS} ${SRC_C}
    ${MEDIAN_SRC_HDRS} ${MEDIAN_SRC_C})
target_add_temporalLogicalRules(${tgt})
target_add_MedianFilter2D_CUDAlib(${tgt})
set_property(TARGET ${tgt} PROPERTY CXX_STANDARD 11)
install(TARGETS ${tgt} RUNTIME DESTINATION test)

if(HAVE_CUDA)
  set(tgt ProcessStack)
  add_executable(${tgt} 
      mainProcessStack.cpp
      ${SRC_CPP} ${SRC_HDRS} ${SRC_C}
      ${MEDIAN_SRC_HDRS} ${MEDIAN_SRC_C})
  target_add_parseconfigfile(${tgt})
  target_add_temporalLogicalRules(${tgt})
  target_add_klb(${tgt})
  target_add_pictools(${tgt})
  target_add_MedianFilter2D_CUDAlib(${tgt})
  set_property(TARGET ${tgt} PROPERTY CXX_STANDARD 11)
  install(TARGETS ${tgt} RUNTIME DESTINATION bin)
endif()

set(tgt ProcessStack_woGPU)
add_executable(${tgt} 
    mainProcessStack.cpp
    ${SRC_CPP} ${SRC_HDRS} ${SRC_C}
    ${MEDIAN_SRC_HDRS} ${MEDIAN_SRC_C})
target_add_parseconfigfile(${tgt})
target_add_temporalLogicalRules_woCUDA(${tgt})
target_add_klb(${tgt})
target_add_pictools(${tgt})
target_add_MedianFilter2D(${tgt})
target_compile_definitions(${tgt} PRIVATE DO_NOT_USE_CUDA)
set_property(TARGET ${tgt} PROPERTY CXX_STANDARD 11)
install(TARGETS ${tgt} RUNTIME DESTINATION bin)

if(WIN32)
    set(tgt ProcessStackBatchMulticore)
    add_executable(${tgt}
        mainProcessStackBatchMulticore.cxx)
    target_add_parseconfigfile(${tgt})
    install(TARGETS ${tgt} RUNTIME DESTINATION bin)
endif()

set(tgt HierarchicalSegmentationPBC)
add_library(${tgt} 
    ${SRC_CPP} ${SRC_HDRS} ${SRC_C}
    ${MEDIAN_SRC_HDRS} ${MEDIAN_SRC_C})
target_add_temporalLogicalRules(${tgt})
target_add_klb(${tgt})
target_add_pictools(${tgt})
target_add_MedianFilter2D_CUDAlib(${tgt})
set_property(TARGET ${tgt} PROPERTY CXX_STANDARD 11)
install(TARGETS ${tgt} ARCHIVE DESTINATION lib)

#organize IDE
source_group("Header Files"  FILES ${SRC_HDRS})
source_group("Header Files\\external\\medianFilter2D"  FILES ${MEDIAN_SRC_HDRS})
source_group("Source Files"  FILES mainTest.cpp mainProcessStack.cpp)
source_group("Source Files"  FILES ${SRC_CPP})
source_group("Source Files"  FILES ${SRC_C})
source_group("Source Files\\external\\medianFilter2D"  FILES ${MEDIAN_SRC_C})


# utility functions

set(HIERARCHICALSEGMENTATION_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

function(target_add_HierarchicalSegmentationPBC tgt)
    target_link_libraries(${tgt} HierarchicalSegmentationPBC)
    target_include_directories(${tgt} PUBLIC ${HIERARCHICALSEGMENTATION_INCLUDE_DIR})
endfunction()
