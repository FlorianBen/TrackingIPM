#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fem.hpp"

int main(int argc, char *argv[]) {

  SpaceCharge::FEMBunch<2> field;

  field.run();

  return 0;
}