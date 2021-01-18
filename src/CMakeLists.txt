# set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)

configure_file("${PROJECT_SOURCE_DIR}/include/SpaceCharge/config.hpp.in"
               "${PROJECT_SOURCE_DIR}/include/SpaceCharge/config.hpp")

set(SPACECHARGE_INC
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/alogger.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/boris.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/bunch.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/config.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/config_reader.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/csv.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/definitions.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/fields.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/nanoflann.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/particle.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/point_cloud.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/run_manager.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/track.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/track_h5.hpp"
    "${PROJECT_SOURCE_DIR}/include/SpaceCharge/units.hpp")

set(SPACECHARGE_SRC "alogger.cpp" "bunch.cpp" "config_reader.cpp" "fields.cpp"
                    "particle.cpp" "run_manager.cpp" "track.cpp")

add_library(SpaceCharge SHARED ${SPACECHARGE_INC} ${SPACECHARGE_SRC})

find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)

if(${USE_MPI})
  find_package(MPI REQUIRED)
  message(
    STATUS
      "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS"
  )
  list(APPEND SPACECHARGE_INC
       "${PROJECT_SOURCE_DIR}/include/SpaceCharge/run_manager_mpi.hpp")
  list(APPEND SPACECHARGE_SRC "run_manager_mpi.cpp")
endif()

if(${USE_HDF5_BLSOC})
  find_package(Blosc REQUIRED)
  list(APPEND SPACECHARGE_INC
       "${PROJECT_SOURCE_DIR}/include/SpaceCharge/blosc_filter.h")
  list(APPEND SPACECHARGE_SRC "blosc_filter.c")
endif()

target_include_directories(
  SpaceCharge
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
  PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)

target_compile_features(SpaceCharge PUBLIC cxx_std_17)
target_link_libraries(SpaceCharge h5cpp GSL::gsl Eigen3::Eigen Boost::boost
                      spdlog::spdlog)

if(${USE_MPI})
  target_link_libraries(SpaceCharge PUBLIC MPI::MPI_CXX)
endif()

if(${USE_HDF5_BLSOC})
  target_link_libraries(SpaceCharge PRIVATE ${BLOSC_LIBRARIES})
endif()

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
