# ##############################################################################
# Sample 001 CMake file
# ##############################################################################
# Author: Florian Benedetti Date: 2021/03/15

set(CLI_SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(TBB REQUIRED)

add_executable(sample_004 ${CLI_SRC})

target_link_libraries(
  sample_004
  PRIVATE SpaceCharge::core
  PRIVATE SpaceCharge::io
  PRIVATE SpaceCharge::field
  PRIVATE SpaceCharge::track
  PRIVATE TBB::tbb)