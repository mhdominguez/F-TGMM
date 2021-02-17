cmake_minimum_required (VERSION 3.5)
project (Tracking_GMM_project)

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)

# C++ Standard
# Notes:
#   17 adds some breaking changes
#       - deletes some previously deprecated interfaces
#       - introduces std::byte which is incompatible with windows headers.
#       - deletes some functions that Eigen used (could maybe update Eigen)
#   14 has experinmental support for <filesystem> which is convenient
#   11 Many features used
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(src/Cmake/git-tag.cmake)
include(src/Cmake/rpath.cmake)
include(src/Cmake/config-klb.cmake)
add_subdirectory(src)
add_subdirectory(doc)

set(CMAKE_INSTALL_DEBUG_LIBRARIES TRUE)
set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
include(InstallRequiredSystemLibraries)
include(CPack) # enable packaging
