add_library(SpaceCharge SHARED "hello.cpp" "bunch.cpp")

configure_file("${PROJECT_SOURCE_DIR}/include/SpaceCharge/version.hpp.in"
               "${PROJECT_BINARY_DIR}/include/SpaceCharge/version.hpp")

find_package(GSL REQUIRED)

target_include_directories(
  SpaceCharge PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                     $<INSTALL_INTERFACE:include>)
target_compile_features(SpaceCharge PUBLIC cxx_std_11)
target_link_libraries(SpaceCharge PUBLIC GSL::gsl)

include(GNUInstallDirs)
install(
  TARGETS SpaceCharge
  EXPORT SpaceChargeTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  #install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/SpaceCharge
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
