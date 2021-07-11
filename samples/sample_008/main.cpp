#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fdm.hpp"

int main(int argc, char *argv[]) {
  SpaceCharge::Logger::Init();

  SpaceCharge::FieldFDM field(600, 400, 0.1, 0.1);

  field.initMatrix();
  field.solve();
  field.save("out.txt");
  return 0;
}