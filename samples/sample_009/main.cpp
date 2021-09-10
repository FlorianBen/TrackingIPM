#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fdm.hpp"
#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/io/field_map_h5.hpp"
#include "SpaceCharge/readout/strips_plane.hpp"

int main(int argc, char *argv[]) {
  using namespace SpaceCharge;
  Logger::Init();

  StripsPlane test(argv[1]);

  quadv<size_t> sizes{0, test.getSizeX(), test.getSizeY(), test.getNbStrips()};
  quadv<double> steps{0, test.getGapX(), test.getGapY(), 1};
  quadv<double> offset{0, 0.0, 0.0, 0.0};

  FieldMap<double> fieldmap(sizes, steps, offset, 0);
  FieldMap<double> potmap(sizes, steps, offset, 0);

  for (auto i = 1; i <= test.getNbStrips(); i++) {
    test.solvePotential(i);
    test.getPotential(potmap, i - 1);
    test.getField(fieldmap, i - 1);
  }

  using namespace hdf5;
  auto file = file::create("out.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();

  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  auto dset1 = root_group.create_dataset(
      "field", datatype::create<FieldMap<double>>(),
      dataspace::create(fieldmap), dcpl, lcpl);
  dset1.write(fieldmap);
  auto dset2 = root_group.create_dataset(
      "pot", datatype::create<FieldMap<double>>(),
      dataspace::create(potmap), dcpl, lcpl);
  dset2.write(potmap);

  return 0;
}