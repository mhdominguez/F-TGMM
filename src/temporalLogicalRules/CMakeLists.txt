add_subdirectory(sparseHungarianAlgorithm)

file(GLOB SRC_HDRS *.h)

file(GLOB XML_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/external/xmlParser/*.h)
file(GLOB XML_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/external/xmlParser/*.cpp)
file(GLOB External_Mylib_HDRS ${MYLIB_INCLUDE_DIR}/*.h)

file(GLOB UtilsCUDA_SRCS knnCUDA/*.cpp)
file(GLOB UtilsCUDA_HDRS knnCUDA/*.h)
file(GLOB UtilsCUDA_CUDA_SRCS knnCUDA/*.cu)

SET(SRC_CPP temporalLogicalRules.cpp GaussianMixtureModel_Redux.cpp nuclei.cpp supervoxel.cpp lineage.cpp lineageHyperTree.cpp trackletCalculation.cpp lineageWindowFeatures.cpp localGeometricDescriptor.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../constants.cpp)

if(HAVE_CUDA)
  CUDA_COMPILE(UtilsCUDA_CUDA_OBJ STATIC ${UtilsCUDA_CUDA_SRCS})
  set(tgt TemporalLogicalRules)
  add_library(${tgt} 
      ${SRC_CPP} ${SRC_HDRS} 
      ${XML_SRCS} ${XML_HDRS} 
      ${UtilsCUDA_SRCS} ${UtilsCUDA_HDRS} 
      ${UtilsCUDA_CUDA_SRCS} 
      ${UtilsExe_SRCS} 
      ${UtilsCUDA_CUDA_OBJ})
  target_add_mylib(${tgt})
  target_add_cuda(${tgt})
  target_add_sparseHungarianAlgorithm(${tgt})
  target_include_directories(${tgt} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
  install(TARGETS ${tgt} ARCHIVE DESTINATION lib)
endif()

set(tgt TemporalLogicalRules_woCUDA)
add_library(${tgt} ${SRC_CPP} ${SRC_HDRS} ${XML_SRCS} ${XML_HDRS} ${UtilsExe_SRCS} )
target_add_mylib(${tgt})
target_add_sparseHungarianAlgorithm(${tgt})
target_include_directories(${tgt} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_definitions(${tgt} PRIVATE DO_NOT_USE_CUDA)
install(TARGETS ${tgt} ARCHIVE DESTINATION lib)

set(TEMPORALLOGICALRULES_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

#organize IDE
source_group("Header Files"  FILES ${SRC_HDRS})
source_group("Header Files\\knnCUDA"  FILES ${UtilsCUDA_HDRS})
source_group("Header Files\\external\\xmlParser"  FILES ${XML_HDRS})
source_group("Source Files\\external\\xmlParser"  FILES ${XML_SRCS})
source_group("Source Files\\knnCUDA"  FILES ${UtilsCUDA_SRCS})
source_group("Source Files\\knnCUDA"  FILES ${UtilsCUDA_CUDA_SRCS})
source_group("Source Files"  FILES ${SRC_CPP})
source_group("Header Files\\mylib" FILES ${External_Mylib_HDRS})


# utility functions for linking to other targets
function(target_add_temporalLogicalRules_woCUDA tgt)
    target_include_directories(${tgt} PUBLIC ${TEMPORALLOGICALRULES_INCLUDE_DIR})
    target_link_libraries(${tgt} TemporalLogicalRules_woCUDA) 
endfunction()

function(target_add_temporalLogicalRules tgt)
    if(HAVE_CUDA)
      target_include_directories(${tgt} PUBLIC ${TEMPORALLOGICALRULES_INCLUDE_DIR})
      target_link_libraries(${tgt} TemporalLogicalRules) 
    else()
      target_add_temporalLogicalRules_woCUDA(${tgt})
    endif()
endfunction()

