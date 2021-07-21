#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fem.hpp"

int main(int argc, char *argv[]) {
  SpaceCharge::Logger::Init();

  SpaceCharge::Particle<double> part("proton", 1, SpaceCharge::cst::mproton,
                                     SpaceCharge::cst::lfactor::beta, 0.5);

  SpaceCharge::FEMBunch<2, double> field(part, 62.5 * SpaceCharge::uni::milli,
          1.0 / (352.0 * SpaceCharge::uni::mega), SpaceCharge::cst::dir::z);

  field.run();

  return 0;
}