cmake_minimum_required(VERSION 3.10)
project(
  SpaceCharge
  VERSION 0.1.0.0
  DESCRIPTION "Calculate Space Charge"
  LANGUAGES CXX C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

option(WITH_BOOST "enable Boost filesystem support" OFF)
set(CONAN "DISABLE" CACHE STRING "conan options AUTO (conan must be in path), MANUAL (expects conanbuildinfo.cmake in build directory) or DISABLE")
add_subdirectory(vendor/h5cpp)

add_subdirectory(src)
add_subdirectory(temp)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

option(BUILD_DOXY_DOC "Build doxygen documentation" ON)
if(BUILD_DOXY_DOC)
  add_subdirectory(docs)
endif()

option(USE_MPI "Enable MPI support." OFF)

option(USE_HDF5_BLSOC "Enable BLOSC compression." OFF)