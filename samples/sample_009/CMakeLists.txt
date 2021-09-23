# ##############################################################################
# Sample 001 CMake file
# ##############################################################################
# Author: Florian Benedetti Date: 2021/03/15

set(CLI_SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_executable(sample_009 ${CLI_SRC})

target_link_libraries(
  sample_009
  PRIVATE SpaceCharge::core
  PRIVATE SpaceCharge::io
  PRIVATE SpaceCharge::field
  PRIVATE SpaceCharge::track
  PRIVATE SpaceCharge::readout)

add_executable(sample_009_bis "main_bis.cpp")

target_link_libraries(
  sample_009_bis
  PRIVATE SpaceCharge::core
  PRIVATE SpaceCharge::io
  PRIVATE SpaceCharge::field
  PRIVATE SpaceCharge::track
  PRIVATE SpaceCharge::readout)
