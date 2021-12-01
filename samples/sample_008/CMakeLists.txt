# ##############################################################################
# Sample 001 CMake file
# ##############################################################################
# Author: Florian Benedetti Date: 2021/03/15

set(CLI_SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_executable(sample_008 ${CLI_SRC})

target_link_libraries(
  sample_008
  PRIVATE SpaceCharge::core
  PRIVATE SpaceCharge::io
  PRIVATE SpaceCharge::field
  PRIVATE SpaceCharge::track
  PRIVATE SpaceCharge::readout)

add_custom_command(
  TARGET sample_008
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/strips_lin.txt
          ${CMAKE_CURRENT_BINARY_DIR}/strips_lin.txt)

add_custom_command(
  TARGET sample_008
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/strips_gaus.txt
          ${CMAKE_CURRENT_BINARY_DIR}/strips_gaus.txt)

add_custom_command(
  TARGET sample_008
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/sonate.txt
          ${CMAKE_CURRENT_BINARY_DIR}/sonate.txt)
