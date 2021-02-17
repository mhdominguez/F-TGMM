file(GLOB GENTLE_BOOST_CPP_SRCS *.cpp)
file(GLOB GENTLE_BOOST_HDRS *.h)

#targets 

add_library (GentleBoost ${GENTLE_BOOST_CPP_SRCS} ${GENTLE_BOOST_HDRS}) 

add_executable (GentleBoost_MainTrain mainTrainSet.cxx ${GENTLE_BOOST_CPP_SRCS} ${GENTLE_BOOST_HDRS}) 

set(tgt GentleBoost_PrecRecallCurve)
add_executable (${tgt}
    mainPrecisionRecallCurve.cxx
    ${GENTLE_BOOST_CPP_SRCS} ${GENTLE_BOOST_HDRS}) 
install(TARGETS ${tgt} RUNTIME DESTINATION bin)


# utilities

set(GENTLEBOOST_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

function(target_add_gentleboost tgt)
    target_include_directories(${tgt} PUBLIC ${GENTLEBOOST_INCLUDE_DIR})
    target_link_libraries(${tgt} GentleBoost)
endfunction()
