file(GLOB CUDA_HDRS *.h)
file(GLOB CUDA_CUDA_SRCS *.cu)

file(GLOB TICTOC_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/external/Nathan/*.h)
file(GLOB TICTOC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/external/Nathan/*.c)

if(HAVE_CUDA)
  CUDA_COMPILE(CUDA_CUDA_OBJ STATIC ${CUDA_CUDA_SRCS} OPTIONS ${CUOPTS})

set(tgt MedianFilter2D)
  add_executable(${tgt} 
      mainTest.cpp
      ${CUDA_SRCS} ${CUDA_HDRS}
      ${CUDA_CUDA_OBJ} ${CUDA_CUDA_SRCS}
      ${TICTOC_HDRS} ${TICTOC_SRCS})
  target_add_cuda(${tgt})
  install(TARGETS ${tgt} RUNTIME DESTINATION test)
endif()

set(tgt MedianFilter2D_CUDAlib)
add_library(${tgt}
    ${CUDA_SRCS} ${CUDA_HDRS}
    ${CUDA_CUDA_OBJ} ${CUDA_CUDA_SRCS}
    ${TICTOC_HDRS} ${TICTOC_SRCS})
target_add_cuda(${tgt})
install(TARGETS ${tgt} ARCHIVE DESTINATION lib)

# utility
set(MEDIANFILTER2DCUDA_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

