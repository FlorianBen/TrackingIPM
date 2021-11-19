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
  /* Initialize logger */
  SpaceCharge::Logger::Init();

  // Create a particle.
  SpaceCharge::Particle<double> part("proton", 1, SpaceCharge::cst::mproton,
                                     SpaceCharge::cst::lfactor::beta, 0.999);

  // Create a bunch
  auto bunch_current = 62.5 * SpaceCharge::uni::milli;
  auto bunch_period = 1.0 / (352.0 * SpaceCharge::uni::mega);
  std::unique_ptr<SpaceCharge::Bunch<double>> bunch(
      new SpaceCharge::GaussianBunch<double>(part, bunch_current, bunch_period,
                                             SpaceCharge::cst::dir::z));

  // Create a generic Field object
  SpaceCharge::FieldSP<double> fields =
      std::make_unique<SpaceCharge::BunchEMField<double>>();
  static_cast<SpaceCharge::BunchEMField<double> *>(fields.get())
      ->usePeriodicity(true);
  static_cast<SpaceCharge::BunchEMField<double> *>(fields.get())
      ->addBunch(std::move(bunch));

  // Create a FieldMap
  size_t nx = 200, ny = 200, nz = 200;
  auto size_x = 50e-3, size_y = 50e-3, size_z = 50e-3;
  SpaceCharge::quadv<size_t> size{0, nx, ny, nz};
  SpaceCharge::quadv<double> step{0, size_x / nx, size_y / ny, size_z / nz};
  SpaceCharge::quadv<double> offset{SpaceCharge::cst::sol * 0e-9, -size_x / 2,
                                    -size_y / 2, -size_z / 2};
  std::unique_ptr<SpaceCharge::FieldMap<double>> fieldmap =
      std::make_unique<SpaceCharge::FieldMapInterpolate<double>>(size, step,
                                                                 offset, 0e-12);
  // Fill FieldMap with our bunch field and compute the EM field at all the
  // given FieldMap coordinates
  static_cast<SpaceCharge::FieldMapInterpolate<double> *>(fieldmap.get())
      ->addField(fields);
  static_cast<SpaceCharge::FieldMapInterpolate<double> *>(fieldmap.get())
      ->computeField();

  using namespace hdf5;
  auto file = file::create("write_fieldmap.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();

  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  // dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  // dcpl.chunk(hdf5::Dimensions{256});

  auto dset = root_group.create_dataset(
      "Efield", datatype::create<SpaceCharge::FieldMap<double>>(),
      dataspace::create(*fieldmap), dcpl, lcpl);
  dset.write(*fieldmap);

  return 0;
}