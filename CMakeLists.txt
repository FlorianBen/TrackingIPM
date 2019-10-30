cmake_minimum_required(VERSION 3.10)
project(
  SpaceCharge
  VERSION 0.1.0.0
  DESCRIPTION "Calculate Space Charge"
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
