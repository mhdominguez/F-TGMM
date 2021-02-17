if(WIN32)
    set(arch "win64")
    set(ext "lib")
    set(pfx "")
else()
    set(arch "linux")
    set(ext "a")
    set(pfx "lib")
endif()

set(MYLIB_LIB_DIR "${CMAKE_CURRENT_LIST_DIR}/../lib/mylib/${arch}/bin")
set(MYLIB_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/../lib/mylib/${arch}/bin")

add_library(mylib STATIC IMPORTED)
set_target_properties(mylib PROPERTIES IMPORTED_LOCATION "${MYLIB_LIB_DIR}/${pfx}mylib.${ext}")

add_library(mytiff STATIC IMPORTED)
set_target_properties(mytiff PROPERTIES IMPORTED_LOCATION "${MYLIB_LIB_DIR}/${pfx}mytiff.${ext}")

function(target_add_mylib tgt)
    target_include_directories(${tgt} PUBLIC ${MYLIB_INCLUDE_DIR})
    target_link_libraries(${tgt} mylib mytiff)
endfunction()