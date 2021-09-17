# ##############################################################################
# CMake file for the readout lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_READOUT
    ${CMAKE_CURRENT_SOURCE_DIR}/ramo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/readout.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/strips_plane.cpp)

set(READOUT_HEADER
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/readout/ramo.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/readout/readout.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/readout/strips_plane.hpp)

install(FILES ${READOUT_HEADER}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/readout)

# Dependencies
find_package(spdlog REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(readout SHARED ${READOUT_HEADER} ${LIB_READOUT})
add_library(SpaceCharge::readout ALIAS readout)

target_link_libraries(
  readout
  PUBLIC core
  PUBLIC track
  PUBLIC spdlog::spdlog
  PUBLIC Eigen3::Eigen)

target_include_directories(
  readout PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
                 $<INSTALL_INTERFACE:include>)

target_compile_features(readout PUBLIC cxx_std_17)
