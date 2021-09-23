#include <iostream>
#include <random>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fdm.hpp"
#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/io/field_map_h5.hpp"
#include "SpaceCharge/readout/ramo.hpp"
#include "SpaceCharge/readout/strips_plane.hpp"

int main(int argc, char *argv[]) {
  using namespace SpaceCharge;
  using namespace hdf5;
  Logger::Init();

  auto map = readMapFromFile<double>(argv[1]);

  auto file = file::create("ramo_cpy.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;

  auto dset_ramo_field =
      root_group.create_dataset("field", datatype::create<FieldMap<double>>(),
                                dataspace::create((*map)), dcpl, lcpl);
  dset_ramo_field.write((*map));

  attribute::Attribute a =
      dset_ramo_field.attributes.create<quadv<double>>("steps");
  quadv<double> steps{0, 1.0, 2.0, 3.0};

  a.write(steps);

  return 0;
}