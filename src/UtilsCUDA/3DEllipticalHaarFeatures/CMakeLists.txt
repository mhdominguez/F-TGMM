add_subdirectory (gentleBoost)

# files

file(GLOB CUDA_HDRS *.h)
file(GLOB CUDA_CUDA_SRCS *.cu)

file(GLOB TICTOC_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/external/Nathan/*.h)
file(GLOB TICTOC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/external/Nathan/*.c)
file(GLOB XML_PARSER_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/external/xmlParser/*.h)
file(GLOB XML_PARSER_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/external/xmlParser/*.cpp)

# targets
if(HAVE_CUDA)
  CUDA_COMPILE(CUDA_CUDA_OBJ STATIC ${CUDA_CUDA_SRCS})


  set(tgt trainCellDivision)
  add_executable(${tgt} 
      mainCellDivisionClassifier.cpp 
      EllipticalHaarFeatures.cpp AnnotationEllipsoid.cpp
      EllipticalHaarFeatures.h AnnotationEllipsoid.h
      ${CUDA_SRCS} ${CUDA_CUDA_OBJ}
      ${CUDA_HDRS} ${CUDA_CUDA_SRCS}
      ${TICTOC_HDRS} ${TICTOC_SRCS}
      ${XML_PARSER_HDRS} ${XML_PARSER_SRCS})
  target_add_cuda(${tgt})
  target_add_mylib(${tgt})
  target_add_gentleboost(${tgt}) 
  install(TARGETS ${tgt} RUNTIME DESTINATION bin)
endif()

set(tgt ellipticalHaarFeatures)
add_library(ellipticalHaarFeatures
    EllipticalHaarFeatures.cpp AnnotationEllipsoid.cpp
    EllipticalHaarFeatures.h AnnotationEllipsoid.h
    ${CUDA_SRCS} ${CUDA_CUDA_OBJ} ${CUDA_HDRS} ${CUDA_CUDA_SRCS}
    ${TICTOC_HDRS} ${TICTOC_SRCS}
    ${XML_PARSER_HDRS} ${XML_PARSER_SRCS})
target_add_cuda(${tgt})
install(TARGETS ${tgt} ARCHIVE DESTINATION lib)

# source groups

source_group("Header Files"  FILES ${CUDA_CUDA_HDRS})
source_group("Header Files\\External\\Nathan"  FILES TICTOC_HDRS)
source_group("Header Files\\External\\xmlParser"  FILES XML_PARSER_HDRS)
source_group("Header Files\\External\\MY_TIFF"  FILES MY_TIFF_HDRS)
source_group("Header Files\\mylib" FILES ${External_Mylib_HDRS})

source_group("Source Files" FILES mainTest.cpp EllipticalHaarFeatures.cpp mainCellDivisionClassifier.cpp AnnotationEllipsoid.cpp ${CUDA_CUDA_SRCS})
source_group("Source Files\\External\\Nathan"  FILES TICTOC_SRCS)
source_group("Source Files\\External\\xmlParser"  FILES XML_PARSER_SRCS)
source_group("Source Files\\External\\MY_TIFF"  FILES MY_TIFF_SRCS)


# utilites

set(ellipticalHaarFeatures_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

function(target_add_ellipticalHaarFeatures tgt)
    target_include_directories(${tgt} PUBLIC ${ellipticalHaarFeatures_INCLUDE_DIR})
    target_link_libraries(${tgt} ellipticalHaarFeatures)
    target_add_gentleboost(${tgt}) 
endfunction()
