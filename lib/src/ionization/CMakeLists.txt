# ##############################################################################
# CMake file for the ionization lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_IONIZATION ${CMAKE_CURRENT_SOURCE_DIR}/bethe_model.cpp)

set(IONIZATION_HEADER
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/ionization/bethe_model.hpp)

install(FILES ${IONIZATION_HEADER}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/ionization)

# Dependencies
find_package(spdlog REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(ionization SHARED ${IONIZATION_HEADER} ${LIB_IONIZATION})
add_library(SpaceCharge::ionization ALIAS ionization)

target_link_libraries(
  ionization
  PUBLIC spdlog::spdlog
  PUBLIC Eigen3::Eigen)

target_include_directories(
  ionization PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
                    $<INSTALL_INTERFACE:include>)

target_compile_features(ionization PUBLIC cxx_std_17)
