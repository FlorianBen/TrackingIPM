# set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)

set(SPACECHARGE_INC
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/bunch.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/particle.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/fields.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/boris.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/csv.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/definitions.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/nanoflann.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/point_cloud.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/track.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/units.hpp")

set(SPACECHARGE_SRC "bunch.cpp" "fields.cpp" "particle.cpp")

add_library(SpaceCharge SHARED ${SPACECHARGE_INC} ${SPACECHARGE_SRC})

configure_file("${PROJECT_SOURCE_DIR}/include/SpaceCharge/version.hpp.in"
               "${PROJECT_BINARY_DIR}/include/SpaceCharge/version.hpp")

find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

target_include_directories(
  SpaceCharge
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
  PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)

target_compile_features(SpaceCharge PUBLIC cxx_std_17)
target_link_libraries(SpaceCharge GSL::gsl Eigen3::Eigen Boost::boost)

include(GNUInstallDirs)
install(
  TARGETS SpaceCharge
  EXPORT SpaceChargeTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  SpaceChargeConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(
  EXPORT SpaceChargeTargets
  FILE SpaceChargeTargets.cmake
  NAMESPACE SpaceCharge::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SpaceCharge)

configure_file(${PROJECT_SOURCE_DIR}/cmake/SpaceChargeConfig.cmake.in
               SpaceChargeConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/SpaceChargeConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/SpaceChargeConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SpaceCharge)
