#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fdm.hpp"
#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/io/field_map_h5.hpp"

int main(int argc, char *argv[]) {
  SpaceCharge::Logger::Init();

  auto size_x = 2000;
  auto size_y = 250;
  auto step = 0.08 / 2000.0;

  SpaceCharge::quadv<size_t> sizes{0, size_x, size_y, 1};
  SpaceCharge::quadv<double> steps{.0, step, step, step};
  SpaceCharge::quadv<double> offset{.0, 0.0, 0.0, 0.0};

  SpaceCharge::FieldMap<double> fieldmap(sizes, steps, offset, .0);
  SpaceCharge::FieldMap<double> potmap(sizes, steps, offset, .0);

  SpaceCharge::FieldFDM field(size_x, size_y, step, step);

  field.readStrips(argv[1]);
  field.initMatrix();
  field.solve();

  field.fillFieldMap(potmap);
  field.fillFieldMap2(fieldmap);

  field.save("out.txt");

  using namespace hdf5;
  auto file = file::create("out.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();

  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  // dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  // dcpl.chunk(hdf5::Dimensions{256});

  auto dset1 = root_group.create_dataset(
      "field", datatype::create<SpaceCharge::FieldMap<double>>(),
      dataspace::create(fieldmap), dcpl, lcpl);
  dset1.write(fieldmap);

  auto dset2 = root_group.create_dataset(
      "pot", datatype::create<SpaceCharge::FieldMap<double>>(),
      dataspace::create(potmap), dcpl, lcpl);
  dset2.write(potmap);

  return 0;
}