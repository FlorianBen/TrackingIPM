# ##############################################################################
# CMake file for tests
# ##############################################################################
# Author: Florian Benedetti Date: 2021/02/02

include(CTest)
enable_testing()
include(GoogleTest)

add_subdirectory("${PROJECT_SOURCE_DIR}/vendor/googletest" "vendor/googletest")

project_add_gtest(test_particle test_particle.cpp SpaceCharge::core "${PROJECT_DIR}")
project_add_gtest(test_bunch test_bunch.cpp SpaceCharge::field "${PROJECT_DIR}")
project_add_gtest(test_field test_field.cpp SpaceCharge::field  "${PROJECT_DIR}")
project_add_gtest(test_logger test_logger.cpp SpaceCharge::core "${PROJECT_DIR}")
