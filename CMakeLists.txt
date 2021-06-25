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

option(USE_MPI "Enable MPI support." OFF)

option(USE_HDF5_BLSOC "Enable BLOSC compression." OFF)

option(USE_DEALII_FEM "Enable support for the deal.II finite element library" OFF)

option(BUILD_DOC "Build doxygen documentation" OFF)

option(BUILD_SAMPLES "Build the samples" ON)

option(BUILD_TESTS "Build the google tests" OFF)

option(ENABLE_COVERAGE "Enable coverage report" OFF)
if(${ENABLE_COVERAGE})
  include(CodeCoverage)
  set(COVERAGE_EXCLUDES
      "${PROJECT_SOURCE_DIR}/vendor/*"
      "${PROJECT_SOURCE_DIR}/tests/*"
      "${PROJECT_SOURCE_DIR}/lib/include/SpaceCharge/nanoflann.hpp"
      "${PROJECT_SOURCE_DIR}/lib/include/SpaceCharge/csv.hpp"
      "/usr/*"
      "/local/home/*/Libs/*")
  append_coverage_compiler_flags()
endif()

option(BUILD_BENCHS "Build the google benchmarks" OFF)

# ##############################################################################
# Add directories according to options
# ##############################################################################

add_subdirectory(lib)

if(BUILD_SAMPLES)
  add_subdirectory(samples)
endif()

if(BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(BUILD_BENCHS)
  add_subdirectory(benchs)
endif()

if(BUILD_DOC)
  add_subdirectory(docs)
endif()

# ##############################################################################
# Coverage report
# ##############################################################################

if(${ENABLE_COVERAGE})
  setup_target_for_coverage_lcov(NAME coverage_spacecharge EXECUTABLE ctest
                                 BASE_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()

# ##############################################################################
# Packaging
# ##############################################################################
set(CPACK_PACKAGE_VENDOR "Florian Benedetti")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "SpaceCharge, a set of tools for particle tracking in an IPM.")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

include(CPack)