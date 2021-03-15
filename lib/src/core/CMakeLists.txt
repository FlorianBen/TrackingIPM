# ##############################################################################
# CMake file for the core lib.
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/24
set(LIB_CORE ${CMAKE_CURRENT_SOURCE_DIR}/alogger.cpp
             ${CMAKE_CURRENT_SOURCE_DIR}/particle.cpp)

set(CORE_HEADER
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/alogger.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/csv.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/definitions.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/nanoflann.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/particle.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/point_cloud.hpp
    ${CMAKE_SOURCE_DIR}/lib/include/SpaceCharge/core/units.hpp)

install(FILES ${CORE_HEADER}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/core)

# Dependencies
find_package(spdlog REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(core SHARED ${CORE_HEADER} ${LIB_CORE})
add_library(SpaceCharge::core ALIAS core)

target_link_libraries(
  core
  PUBLIC spdlog::spdlog
  PUBLIC Eigen3::Eigen)

target_include_directories(
  core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include/>
              $<INSTALL_INTERFACE:include>)

target_compile_features(core PUBLIC cxx_std_17)
