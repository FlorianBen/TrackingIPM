# ##############################################################################
# CMake file for libs
# ##############################################################################
# Author: Florian Benedetti Date: 2020/03/12
add_subdirectory(src/core)
add_subdirectory(src/io)
add_subdirectory(src/field)
add_subdirectory(src/track)

configure_file("${PROJECT_SOURCE_DIR}/lib/include/SpaceCharge/config.hpp.in"
               "${PROJECT_SOURCE_DIR}/lib/include/SpaceCharge/config.hpp")
install(FILES "${PROJECT_SOURCE_DIR}/lib/include/SpaceCharge/config.hpp"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpaceCharge/)

include(GNUInstallDirs)
install(
  TARGETS core io field track
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
