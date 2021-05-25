#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/field/fields.hpp"
#include "SpaceCharge/io/field_map_h5.hpp"

int main(int argc, char *argv[]) {
  /* code */
  SpaceCharge::Logger::Init();

  SpaceCharge::Particle<double> part("proton", 1, SpaceCharge::cst::mproton,
                                     SpaceCharge::cst::lfactor::beta, 0.5);

  std::unique_ptr<SpaceCharge::Bunch<double>> bunch2(
      new SpaceCharge::GaussianBunch<double>(
          part, 62.5 * SpaceCharge::uni::milli,
          1.0 / (352.0 * SpaceCharge::uni::mega), SpaceCharge::cst::dir::z));

  SpaceCharge::FieldSP<double> fields =
      std::make_unique<SpaceCharge::BunchEMField<double>>();
  static_cast<SpaceCharge::BunchEMField<double> *>(fields.get())
      ->usePeriodicity(true);
  static_cast<SpaceCharge::BunchEMField<double> *>(fields.get())
      ->addBunch(std::move(bunch2));

  auto nx = 400;
  auto ny = 400;
  auto nz = 1;

  using namespace hdf5;
  auto file = file::create("write_fieldmap.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  auto data_dir = root_group.create_group("data");

  SpaceCharge::quadv<size_t> size{0, nx, ny, nz};
  SpaceCharge::quadv<double> step{0, 100.0e-3 / nx, 100.0e-3 / ny,
                                  100.0e-3 / nz};

  for (auto i = 1; i < 2; i++) {
    SpaceCharge::quadv<double> offset{SpaceCharge::cst::sol * 1e-9 * i, -50e-3,
                                      -50e-3, -0e-3};

    SpaceCharge::FieldMap<double> fieldmap(size, step, offset, 3e-12);
    fieldmap.addField(fields);
    fieldmap.computeField();

    auto dset = data_dir.create_dataset(
        "sequence_" + std::to_string(i),
        datatype::create<SpaceCharge::FieldMap<double>>(),
        dataspace::create(fieldmap), dcpl, lcpl);
    dset.write(fieldmap);
  }

  return 0;
}