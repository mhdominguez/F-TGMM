# files

file(GLOB LIB_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB LIB_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB TICTOC_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/../external/Nathan/*.h)
file(GLOB TICTOC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../external/Nathan/*.c)
file(GLOB XML_PARSER_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/../external/xmlParser2/*.h)
file(GLOB XML_PARSER_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../external/xmlParser2/*.cpp)

# targets
if(HAVE_CUDA)
  set(tgt trainCellDivisionWithTemporalWindow)
  add_executable(${tgt} mainCellDivisionClassifierWithTemporalWindow.cxx
      ${LIB_HDRS} ${LIB_SRCS}
      ${TICTOC_HDRS} ${TICTOC_SRCS}
      ${XML_PARSER_HDRS} ${XML_PARSER_SRCS})
  target_add_mylib(${tgt})
  target_add_temporalLogicalRules(${tgt})
  target_add_ellipticalHaarFeatures(${tgt})
  target_add_pictools(${tgt})
  target_add_klb(${tgt})
  install(TARGETS ${tgt} RUNTIME DESTINATION bin)
endif()

set(tgt trainCellDivisionWithTemporalWindowLib)
add_library(${tgt}
    ${LIB_HDRS} ${LIB_SRCS}
    ${TICTOC_HDRS} ${TICTOC_SRCS}
    ${XML_PARSER_HDRS} ${XML_PARSER_SRCS})
target_add_mylib(${tgt})
target_add_temporalLogicalRules(${tgt})
target_add_ellipticalHaarFeatures(${tgt})
target_add_pictools(${tgt})
target_add_klb(${tgt})
install(TARGETS ${tgt} ARCHIVE DESTINATION lib)

# source groups

source_group("Header Files\\External\\Nathan"  FILES TICTOC_HDRS)
source_group("Header Files\\External\\xmlParser"  FILES XML_PARSER_HDRS)
source_group("Header Files\\External\\MY_TIFF"  FILES MY_TIFF_HDRS)

source_group("Source Files" FILES mainCellDivisionClassifierWithTemporalWindow.cxx)
source_group("Source Files\\External\\Nathan"  FILES TICTOC_SRCS)
source_group("Source Files\\External\\xmlParser"  FILES XML_PARSER_SRCS)
source_group("Source Files\\External\\MY_TIFF"  FILES MY_TIFF_SRCS)

#utilities

set(trainCellDivisionWithTemporalWindowLib_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

function(target_add_trainCellDivisionWithTemporalWindowLib tgt)
    target_include_directories(${tgt} PUBLIC ${trainCellDivisionWithTemporalWindowLib_INCLUDE_DIR})
    target_link_libraries(${tgt} trainCellDivisionWithTemporalWindowLib)
endfunction()
