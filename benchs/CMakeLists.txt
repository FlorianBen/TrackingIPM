# ##############################################################################
# CMake file for benchmarks
# ##############################################################################
# Author: Florian Benedetti Date: 2021/06/13

add_subdirectory("${PROJECT_SOURCE_DIR}/vendor/benchmark" "vendor/benchmark")

add_executable(bench_quadv bench_quadv.cpp)
target_link_libraries(bench_quadv SpaceCharge::core benchmark::benchmark)