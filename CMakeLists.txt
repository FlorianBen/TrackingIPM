cmake_minimum_required(VERSION 3.10)
project(
  SpaceCharge
  VERSION 0.1.0.0
  DESCRIPTION "Calculate Space Charge"
  LANGUAGES CXX C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

add_subdirectory(src)

add_subdirectory(temp)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

option(BUILD_DOXY_DOC "build doxygen documentation" OFF)
if(BUILD_DOXY_DOC)
  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(WARNING "Doxygen not found, unable to generate documentation")
  else()
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")
    doxygen_add_docs(libSpaceChargeDoc ${CMAKE_SOURCE_DIR}/include)
  endif()
endif()
