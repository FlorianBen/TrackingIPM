# ##############################################################################
# CMake file for the io lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_IO ${CMAKE_CURRENT_SOURCE_DIR}/blosc_filter.c)

set(HEADER_IO
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/io/blosc_filter.h
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/io/field_map_h5.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/io/quadv_h5.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/io/track_h5.hpp)

install(FILES ${HEADER_IO}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/io)

# Dependencies
add_subdirectory(${CMAKE_SOURCE_DIR}/vendor/h5cpp build_h5cpp)

add_library(io SHARED ${HEADER_IO} ${LIB_IO})
add_library(SpaceCharge::io ALIAS io)

target_link_libraries(
  io
  PUBLIC core
  PUBLIC h5cpp)

target_include_directories(
  io PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
            $<INSTALL_INTERFACE:include>)

target_compile_features(io PUBLIC cxx_std_17)
