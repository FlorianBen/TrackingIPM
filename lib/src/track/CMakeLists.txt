# ##############################################################################
# CMake file for the io lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_TRACK ${CMAKE_CURRENT_SOURCE_DIR}/track.cpp)

set(HEADER_TRACK ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/track/track.hpp)

install(FILES ${HEADER_TRACK}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/track/)

# Dependencies
find_package(Boost REQUIRED)

add_library(track SHARED ${HEADER_TRACK} ${LIB_TRACK})
add_library(SpaceCharge::track ALIAS track)

target_link_libraries(
  track
  PUBLIC core
  PUBLIC Boost::boost)

target_include_directories(
  track PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
               $<INSTALL_INTERFACE:include>)

target_compile_features(io PUBLIC cxx_std_17)
