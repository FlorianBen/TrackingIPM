# ##############################################################################
# Main CMake file
# ##############################################################################
# Author: Florian Benedetti Date: 2021/02/02

cmake_minimum_required(VERSION 3.10)
project(
  SpaceCharge
  VERSION 0.1.0.0
  DESCRIPTION "Calculate Space Charge"
  LANGUAGES CXX C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(SpaceChargeUtils)

# ##############################################################################
# Options and configuration
# ##############################################################################

option(WITH_BOOST "enable Boost filesystem support" OFF)
set(CONAN
    "DISABLE"
    CACHE
      STRING
      "conan options AUTO (conan must be in path), MANUAL (expects conanbuildinfo.cmake in build directory) or DISABLE"
)
add_subdirectory(vendor/h5cpp)

option(USE_MPI "Enable MPI support." OFF)

option(USE_HDF5_BLSOC "Enable BLOSC compression." OFF)

option(BUILD_DOC "Build doxygen documentation" OFF)

option(BUILD_TESTS "Build the google tests" OFF)


# ##############################################################################
# Add directories according to options
# ##############################################################################

add_subdirectory(src)
add_subdirectory(temp)

if(BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(BUILD_DOC)
  add_subdirectory(docs)
endif()

option(USE_MPI "Enable MPI support." OFF)

option(USE_HDF5_BLSOC "Enable BLOSC compression." OFF)