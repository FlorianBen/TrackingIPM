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
  typedef FieldMap<double> fmap;
  typedef std::unique_ptr<fmap> upfmap;
  typedef std::unique_ptr<Track<double>> uptrack;

  Logger::Init();

  using namespace hdf5;
  auto file = file::create("ramo.h5", file::AccessFlags::TRUNCATE);
  auto root_group = file.root();
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;

  SpaceCharge::Particle<double> part_track(
      "electron", 1, 2.0 * SpaceCharge::cst::mproton,
      SpaceCharge::cst::lfactor::beta, 0.5);

  std::default_random_engine g1(0);                  // use always same 0 seed.
  std::normal_distribution<double> pos_x(0.0, 3e-3); // pos x distribution
  std::normal_distribution<double> pos_y(0.0, 3e-3); // pos y distribution
  std::uniform_real_distribution<double> pos_z(-0.012, 0.012);
  std::normal_distribution<double> pos_t(-2e-12, 2e-12);
  auto nb_part = 1;
  tbb::concurrent_vector<SpaceCharge::quadv<double>> pos(nb_part);
  for (auto &p : pos) {
    p(0) = 0.0;
    p(1) = 0.0; // p(1) = -0.020; // pos_x(g1);
    p(2) = 0.0; // pos_y(g1);
    p(3) = 0.0; // pos_z(g1);
  }
  SpaceCharge::quadv<double> v0{0.0, 0.0, 0.0e5, 0.0};

  SpaceCharge::FieldSPS<double> fields =
      std::make_shared<SpaceCharge::EMFieldsManager<double>>();
  SpaceCharge::quadv<double> Ecst{0.0, 0.0, 3.0e5, 0.0};
  SpaceCharge::quadv<double> Bcst{0.0, 0.0, 0.0, 0.0};
  SpaceCharge::state_type2<double> EMcst{Ecst, Bcst};
  SpaceCharge::FieldSP<double> Fep =
      std::make_unique<SpaceCharge::ConstantEMField<double>>(EMcst);
  static_cast<SpaceCharge::EMFieldsManager<double> *>(fields.get())
      ->addField(Fep);

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

  TrackSPS<double> tr = std::make_shared<Track<double>>(tracks.back());

  StripsPlane test(argv[1]);
  quadv<size_t> sizes{0, test.getSizeX(), test.getSizeY(), test.getNbStrips()};
  quadv<double> steps{0, test.getGapX(), test.getGapY(), 1};
  quadv<double> offset{0, test.getGapX() * test.getSizeX() / 2.0, 0.0, 0.0};

  upfmap temp = std::make_unique<fmap>(sizes, steps, offset, 0);
  upfmap ramopot = std::make_unique<fmap>(sizes, steps, offset, 0);

  FieldMapSPS<double> ramofield = std::move(temp);

  for (auto i = 1; i <= test.getNbStrips(); i++) {
    test.solvePotential(i);
    test.getPotential(*ramopot, i - 1);
    test.getField(*ramofield, i - 1);
  }

  auto dset_ramo_field =
      root_group.create_dataset("field", datatype::create<FieldMap<double>>(),
                                dataspace::create((*ramofield)), dcpl, lcpl);

  attribute::Attribute a = dset_ramo_field.attributes.create<quadv<double>>("note");
  a.write(steps);

  dset_ramo_field.write((*ramofield));

  auto dset_ramo_pot =
      root_group.create_dataset("pot", datatype::create<FieldMap<double>>(),
                                dataspace::create((*ramopot)), dcpl, lcpl);
  dset_ramo_pot.write((*ramopot));

  RamoComputation<double> ramo(ramofield, tr);

  auto dset_ramo_coor = root_group.create_dataset(
      "pos", datatype::create<std::vector<SpaceCharge::quadv<int>>>(),
      dataspace::create(ramo.getCoordinates()), dcpl, lcpl);
  dset_ramo_coor.write(ramo.getCoordinates());

  auto dset_ramo_traj = root_group.create_dataset(
      "traj", datatype::create<std::vector<SpaceCharge::quadv<double>>>(),
      dataspace::create(ramo.getTrajectory()), dcpl, lcpl);
  dset_ramo_traj.write(ramo.getTrajectory());

  auto current_group = root_group.create_group("current");
  auto i = 0;
  for (const auto &cv : ramo.getCurrent()) {
    auto dset_ramo_current =
        current_group.create_dataset("current_" + std::to_string(i++),
                                     datatype::create<std::vector<double>>(),
                                     dataspace::create(cv), dcpl, lcpl);
    dset_ramo_current.write(cv);
  }

  return 0;
}