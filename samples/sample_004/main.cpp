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

  size_t nx = 1;
  size_t ny = 1;
  size_t nz = 400;

  SpaceCharge::quadv<size_t> size{0, nx, ny, nz};
  SpaceCharge::quadv<double> step{0, 100.0e-3 / nx, 100.0e-3 / ny,
                                  100.0e-3 / nz};
  SpaceCharge::quadv<double> offset{SpaceCharge::cst::sol * 0e-9, -0e-3, -0e-3,
                                    -50e-3};

  std::unique_ptr<SpaceCharge::FieldMap<double>> fieldmap =
      std::make_unique<SpaceCharge::FieldMapInterpolate<double>>(size, step,
                                                                 offset, 3e-12);
  static_cast<SpaceCharge::FieldMapInterpolate<double> *>(fieldmap.get())->addField(fields);
  static_cast<SpaceCharge::FieldMapInterpolate<double> *>(fieldmap.get())->computeField();

  using namespace hdf5;
  auto file = file::create("write_fieldmap.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();

  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  // dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  // dcpl.chunk(hdf5::Dimensions{256});

  auto dset = root_group.create_dataset(
      "sequence", datatype::create<SpaceCharge::FieldMap<double>>(),
      dataspace::create(*fieldmap), dcpl, lcpl);
  dset.write(*fieldmap);

  return 0;
}