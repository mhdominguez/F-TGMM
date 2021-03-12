cmake_minimum_required (VERSION 2.8)
project (mylib)

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

add_subdirectory(src)

include(InstallRequiredSystemLibraries)
include(CPack) # enable packaging
