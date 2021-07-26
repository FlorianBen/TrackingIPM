#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fdm.hpp"

int main(int argc, char *argv[]) {
  SpaceCharge::Logger::Init();

  SpaceCharge::FieldFDM field(2000, 250, 0.08/2000, 0.01/500);

  field.readStrips(argv[1]);
  field.initMatrix();
  field.solve();
  field.save("out.txt");
  return 0;
}