#include <boost/numeric/odeint.hpp>
#include <h5cpp/hdf5.hpp>
#include <iostream>
#include <random>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/core/point_cloud.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/fields.hpp"
#include "SpaceCharge/io/track_h5.hpp"
#include "SpaceCharge/track/track.hpp"

int main(int argc, char *argv[]) {
  SpaceCharge::Logger::Init();

  // Create two particles. Each can be easily modified by changing the
  // constructor arguments.
  SpaceCharge::Particle<double> part_bunch(
      "proton", 1, 1.0 * SpaceCharge::cst::mproton,
      SpaceCharge::cst::lfactor::beta, 0.5);
  SpaceCharge::Particle<double> part_track(
      "proton", 1, 1.0 * SpaceCharge::cst::mproton,
      SpaceCharge::cst::lfactor::beta, 0.5);

  // Setup the random number generator.
  std::default_random_engine g1(0);                  // use always same 0 seed.
  std::normal_distribution<double> pos_x(0.0, 3e-3); // pos x distribution
  std::normal_distribution<double> pos_y(0.0, 3e-3); // pos y distribution
  std::uniform_real_distribution<double> pos_z(-0.012, 0.012);
  std::normal_distribution<double> pos_t(-2e-12, 2e-12);

  SpaceCharge::quadv<double> v0{0.0, 0.0, 0.0, 0.0};
  auto t0 = 0.0;

  SpaceCharge::FieldSPS<double> fields =
      std::make_shared<SpaceCharge::EMFieldsManager<double>>();

  SpaceCharge::quadv<double> offset;
  offset << 0.0, 0.0, 0.0, 0.206;

  SpaceCharge::FieldSP<double> Fep =
      std::make_unique<SpaceCharge::CSVFileEMField<double>>();
  static_cast<SpaceCharge::CSVFileEMField<double> *>(Fep.get())->loadEfield(
      "DataRaw.csv", offset);

  // SpaceCharge::quadv<double> Ecst{0.0, 3.0e5, 0.0, 0.0};
  // SpaceCharge::quadv<double> Bcst{0.0, 0.0e5, 0.0, 0.0};
  // SpaceCharge::state_type2<double> EMcst{Ecst, Bcst};
  // SpaceCharge::FieldSP<double> Fep =
  //    std::make_unique<SpaceCharge::ConstantEMField<double>>(EMcst);

  static_cast<SpaceCharge::EMFieldsManager<double> *>(fields.get())
      ->addField(Fep);

  //   std::unique_ptr<SpaceCharge::Bunch<double>> bunch2(
  //       new SpaceCharge::GaussianBunch<double>(
  //           part_bunch, 62.5 * SpaceCharge::uni::milli,
  //           1.0 / (352.0 * SpaceCharge::uni::mega),
  //           SpaceCharge::cst::dir::z));
  //   SpaceCharge::FieldSP<double> field_bunch =
  //       std::make_unique<SpaceCharge::BunchEMField<double>>();
  //   static_cast<SpaceCharge::BunchEMField<double> *>(field_bunch.get())
  //       ->usePeriodicity(true);
  //   static_cast<SpaceCharge::BunchEMField<double> *>(field_bunch.get())
  //       ->addBunch(std::move(bunch2));
  //   static_cast<SpaceCharge::EMFieldsManager<double> *>(fields.get())
  //       ->addField(field_bunch);

  auto nb_part = 4000;
  tbb::concurrent_vector<SpaceCharge::quadv<double>> pos(nb_part);
  for (auto &p : pos) {
    p(0) = 0.0;
    p(1) = pos_x(g1);
    p(2) = pos_y(g1);
    p(3) = pos_z(g1);
  }

  tbb::concurrent_vector<SpaceCharge::Track<double>> tracks;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, nb_part),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (auto i = r.begin(); i != r.end(); i++) {
                        if ((i % 1000) == 0) {
                          SC_INFO("Done {}/{}", i, nb_part);
                        }
                        SpaceCharge::Track<double> track(part_track, pos[i], v0,
                                                         fields);
                        track.track();
                        tracks.push_back(track);
                      }
                    });

  auto i = 0;

  SC_INFO("Tracking done");
  SC_INFO("Saving");

  std::string file = "tracking.h5";
  hdf5::file::File f =
      hdf5::file::create(file, hdf5::file::AccessFlags::TRUNCATE);

  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  // dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  // dcpl.chunk(hdf5::Dimensions{256});

  auto root_group = f.root();

  using namespace hdf5;

  for (auto &track : tracks) {
    node::Group group =
        root_group.create_group("track_" + std::to_string(i++), lcpl);
    auto dset = group.create_dataset(
        "pos", datatype::create<std::vector<SpaceCharge::quadv<double>>>(),
        dataspace::create(track.getPosVector()), dcpl, lcpl);
    dset.write(track.getPosVector());
    auto dset_speed = group.create_dataset(
        "speed", datatype::create<std::vector<SpaceCharge::quadv<double>>>(),
        dataspace::create(track.getSpeedVector()), dcpl, lcpl);
    dset_speed.write(track.getSpeedVector());
  }

  SC_INFO("Saving done");

  return 0;
}