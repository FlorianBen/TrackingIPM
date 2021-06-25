# ##############################################################################
# CMake file for the io lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_FIELD
    ${CMAKE_CURRENT_SOURCE_DIR}/bunch.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/field_map.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/field_fdm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fields.cpp)

set(HEADER_FIELD
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/bunch.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_map.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_comsol.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/field/field_map.hpp)

install(FILES ${HEADER_FIELD}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/field/)

# Dependencies
find_package(TBB REQUIRED)
find_package(GSL REQUIRED)

add_library(field SHARED ${HEADER_FIELD} ${LIB_FIELD})
add_library(SpaceCharge::field ALIAS field)

target_link_libraries(
  field
  PUBLIC core
  PUBLIC TBB::tbb
  PUBLIC GSL::gsl)

target_include_directories(
  field PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
               $<INSTALL_INTERFACE:include>)

target_compile_features(io PUBLIC cxx_std_17)
