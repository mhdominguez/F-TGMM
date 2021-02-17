#Cuda
#You might have to deactive any device debugging symbols (-G) to avoid a kernel launch failure in thurst::sort_by_key
#SET(CUDA_NVCC_FLAGS "--cudart static;-gencode arch=compute_30,code=sm_30;-gencode arch=compute_32,code=sm_32;-gencode arch=compute_35,code=sm_35;-gencode arch=compute_50,code=sm_50;-gencode arch=compute_50,code=compute_50" CACHE STRING "Semi-colon delimit multiple arguments")
SET(CUDA_NVCC_FLAGS "--cudart static;-gencode arch=compute_52,code=sm_52;-gencode arch=compute_60,code=sm_60;-gencode arch=compute_61,code=sm_61;-gencode arch=compute_61,code=compute_61" CACHE STRING "Semi-colon delimit multiple arguments")
SET(CUDA_NVCC_FLAGS_DEBUG -g; -G;-O0;-DTHRUST_DEBUG CACHE STRING "Semi-colon delimit multiple arguments") #set before FIND_PACKAGE(CUDA) in order to avoid FORCE to show them in GUI. So user can modify them
SET(CUDA_NVCC_FLAGS_RELEASE -O3 CACHE STRING "Semi-colon delimit multiple arguments")
SET(CUDA_NVCC_FLAGS_RELWITHDEBINFO -O3 CACHE STRING "Semi-colon delimit multiple arguments")
FIND_PACKAGE(CUDA) # To select the cuda version, set CUDA_BIN_PATH appropriately in your environment. Or set CUDA_TOOLKIT_ROOT_DIR after configuration.
SET(CUDA_VERBOSE_BUILD ON FORCE)

file(GLOB CUDA_SHARED_LIBS "${CUDA_TOOLKIT_ROOT_DIR}/bin/*.dll")
install(FILES ${CUDA_SHARED_LIBS} DESTINATION bin)

set(HAVE_CUDA ${CUDA_FOUND})

function(target_add_cuda tgt)    
    if(HAVE_CUDA) 
      target_link_libraries(${tgt} 
          ${CUDA_cusparse_LIBRARY}
          ${CUDA_CUDART_LIBRARY})
      CUDA_ADD_CUBLAS_TO_TARGET(${tgt})
      CUDA_ADD_CUFFT_TO_TARGET(${tgt})
      target_include_directories(${tgt} PUBLIC ${CUDA_INCLUDE_DIRS})
    else()
      target_compile_definitions(${tgt} PUBLIC DO_NOT_USE_CUDA)
    endif()
endfunction()
