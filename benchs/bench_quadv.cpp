#include "SpaceCharge/core/point_cloud.hpp"
#include <benchmark/benchmark.h>

static void BM_QuadvCreation(benchmark::State &state) {
  for (auto _ : state)
    SpaceCharge::quadv<double> tquadv;
}
// Register the function as a benchmark
BENCHMARK(BM_QuadvCreation);

static void BM_QuadvAddition(benchmark::State &state) {
  SpaceCharge::quadv<double> tquadv1{0.0, 1.0, 2.0, 3.0};
  SpaceCharge::quadv<double> tquadv2{1.0, 2.0, 3.0, 4.0};
  for (auto _ : state)
    tquadv1 + tquadv2;
}
// Register the function as a benchmark
BENCHMARK(BM_QuadvAddition);

static void BM_QuadvScalarProd(benchmark::State &state) {
  SpaceCharge::quadv<double> tquadv1{0.0, 1.0, 2.0, 3.0};
  SpaceCharge::quadv<double> tquadv2{1.0, 2.0, 3.0, 4.0};
  for (auto _ : state)
    SpaceCharge::scalar_prod(tquadv1,tquadv2);
}
// Register the function as a benchmark
BENCHMARK(BM_QuadvScalarProd);

static void BM_QuadvVectorProduct(benchmark::State &state) {
  SpaceCharge::quadv<double> tquadv1{0.0, 1.0, 2.0, 3.0};
  SpaceCharge::quadv<double> tquadv2{1.0, 2.0, 3.0, 4.0};
  for (auto _ : state)
    SpaceCharge::vect_prod(tquadv1,tquadv2);
}
// Register the function as a benchmark
BENCHMARK(BM_QuadvVectorProduct);

BENCHMARK_MAIN();