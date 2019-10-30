add_library(SpaceCharge SHARED "bunch.cpp" "particle.cpp")

configure_file("${PROJECT_SOURCE_DIR}/include/SpaceCharge/version.hpp.in"
               "${PROJECT_BINARY_DIR}/include/SpaceCharge/version.hpp")

find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

target_include_directories(
  SpaceCharge PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                     $<INSTALL_INTERFACE:include>)
target_compile_features(SpaceCharge PUBLIC cxx_std_14)
target_link_libraries(SpaceCharge PUBLIC GSL::gsl)
target_link_libraries(SpaceCharge PUBLIC Boost::boost)

include(GNUInstallDirs)
install(
  TARGETS SpaceCharge
  EXPORT SpaceChargeTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
