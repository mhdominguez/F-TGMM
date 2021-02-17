include_directories(${CMAKE_CURRENT_LIST_DIR})

include(Cmake/config-pictools.cmake)
include(Cmake/config-cuda-support.cmake)
include(Cmake/config-mylib.cmake)
include(Cmake/config-parseconfigfile.cmake)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory (temporalLogicalRules)
add_subdirectory (nucleiChSvWshedPBC)
add_subdirectory (UtilsCUDA/3DEllipticalHaarFeatures)
add_subdirectory (backgroundDetection)
add_subdirectory (cellDivisionDiscriminationWithTemporalWindow)

### TGMM Executable 

file(GLOB UtilsExe_SRCS UtilsExe/*.cpp)
file(GLOB UtilsExe_HDRS UtilsExe/*.h)

file(GLOB Utils_SRCS Utils/*.cpp)
file(GLOB UtilsC_SRCS Utils/*.c)
file(GLOB Utils_HDRS Utils/*.h)

if(HAVE_CUDA)
  file(GLOB UtilsCUDA_CUDA_SRCS UtilsCUDA/*.cu)
  file(GLOB UtilsCUDA_SRCS UtilsCUDA/*.cpp)
  file(GLOB UtilsCUDA_HDRS UtilsCUDA/*.h)
  CUDA_COMPILE(UtilsCUDA_CUDA_OBJ STATIC ${UtilsCUDA_CUDA_SRCS} constants.cpp)
endif()


file(GLOB_RECURSE ExternalC_SRCS external/*.c)
file(GLOB_RECURSE ExternalCXX_SRCS external/*.cpp)
file(GLOB_RECURSE External_HDRS external/*.h)

file(GLOB ConfigFiles filesToCopyToExeFolder/*.*)
install(FILES ${ConfigFiles} DESTINATION bin)

set(MAIN_CLASSES_SRCS GaussianMixtureModel.cpp variationalInference.cpp kdtree.cpp responsibilities.cpp cellDivision.cpp supportFunctionsMain.cpp  backgroundDetectionInterface.cpp constants.cpp)
file(GLOB MAIN_CLASSES_HDRS *.h)

if(HAVE_CUDA)
  add_executable(TGMM mainTrackingGaussianMixtureModel.cpp 
      ${MAIN_CLASSES_SRCS} 
      ${MAIN_CLASSES_HDRS} 
      ${UtilsCUDA_SRCS} 
      ${UtilsCUDA_CUDA_OBJ} 
      ${UtilsCUDA_HDRS}
      ${UtilsCUDA_CUDA_SRCS}
      ${UtilsExe_SRCS}
      ${UtilsExe_HDRS}
      ${Utils_SRCS}
      ${UtilsC_SRCS}
      ${Utils_HDRS} 
      ${OpticalFlow_SRCS} 
      ${OpticalFlow_HDRS} 
      ${UtilsMatlabSharedLibrary_SRCS} 
      ${UtilsMatlabSharedLibrary_HDRS} 
      ${ExternalC_SRCS} 
      ${ExternalCXX_SRCS} 
      ${External_HDRS})
  target_add_cuda(TGMM)
  target_add_mylib(TGMM)
  target_add_pictools(TGMM)
  target_add_klb(TGMM)
  target_add_temporalLogicalRules(TGMM)
  target_add_BackgroundDetectionLib(TGMM)
  target_add_ellipticalHaarFeatures(TGMM)
  target_add_HierarchicalSegmentationPBC(TGMM)
  target_add_trainCellDivisionWithTemporalWindowLib(TGMM)
  set_property(TARGET TGMM PROPERTY CXX_STANDARD 11)
  install(TARGETS TGMM RUNTIME DESTINATION bin)
endif()

# 
# Source groups for IDE
# 

#We need each external folder separately
file(GLOB External_gsl_SRCS external/gsl/*.c)
file(GLOB External_Nathan_SRCS external/Nathan/*.c)
file(GLOB External_xmlParser_SRCS external/xmlParser2/*.cpp)

file(GLOB External_Mylib_HDRS mylib/*.h)
file(GLOB_RECURSE External_Eigen_HDRS external/Eigen/*)

file(GLOB External_gsl_HDRS external/gsl/*.h)
file(GLOB External_Nathan_HDRS external/Nathan/*.h)
file(GLOB External_xmlParser_HDRS external/xmlParser2/*.h)


source_group("Header Files"  FILES ${MAIN_CLASSES_HDRS})
source_group("Header Files\\UtilsCUDA"  FILES ${UtilsCUDA_HDRS})
source_group("Header Files\\Utils"  FILES ${Utils_HDRS})
source_group("Header Files\\External\\gsl" FILES ${External_gsl_HDRS})
source_group("Header Files\\External\\Nathan" FILES ${External_Nathan_HDRS})
source_group("Header Files\\External\\xmlParser" FILES ${External_xmlParses_HDRS})
source_group("Header Files\\mylib" FILES ${External_Mylib_HDRS})
source_group("Header Files\\External\\Eigen" FILES ${External_Eigen_HDRS})

source_group("Source Files" FILES mainTrackingGaussianMixtureModel.cpp)
source_group("Source Files\\UtilsCUDA"  FILES ${UtilsCUDA_SRCS})
source_group("Source Files\\UtilsCUDA"  FILES ${UtilsCUDA_CUDA_SRCS})
source_group("Source Files\\Utils"  FILES ${Utils_SRCS})
source_group("Source Files\\Utils"  FILES ${UtilsC_SRCS})
source_group("Source Files\\External\\gsl" FILES ${External_gsl_SRCS})
source_group("Source Files\\External\\Nathan" FILES ${External_Nathan_SRCS})
source_group("Source Files\\External\\xmlParser" FILES ${External_xmlParses_SRCS})
