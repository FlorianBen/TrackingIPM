#include <atomic>
#include <iostream>
#include <random>
#include <thread>
#include <utility>

#include <h5cpp/hdf5.hpp>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/field_fdm.hpp"
#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/io/field_map_h5.hpp"
#include "SpaceCharge/readout/ramo.hpp"
#include "SpaceCharge/readout/strips_plane.hpp"

using ramoSP = std::shared_ptr<SpaceCharge::RamoComputation<double>>;
tbb::concurrent_bounded_queue<ramoSP> ramos;
std::atomic_bool prod = true;

void save() {
  using namespace hdf5;
  auto file_out = file::create("current.h5", file::AccessFlags::TRUNCATE);
  auto root_group_out = file_out.root();
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  ramoSP res;
  auto i = 0;
  SC_INFO("Saving loop");
  while (!ramos.empty() || prod) {
    auto new_data = ramos.try_pop(res);
    if (new_data) {
      std::string ind_part_sufix =
          std::string(5 - std::to_string(i).length(), '0') + std::to_string(i);
      i++;
      auto part_group = root_group_out.create_group("part_" + ind_part_sufix);
      auto dset_ramo_coor = part_group.create_dataset(
          "pos", datatype::create<std::vector<SpaceCharge::quadv<int>>>(),
          dataspace::create(res->getCoordinates()), dcpl, lcpl);
      dset_ramo_coor.write(res->getCoordinates());

      auto dset_ramo_traj = part_group.create_dataset(
          "traj", datatype::create<std::vector<SpaceCharge::quadv<double>>>(),
          dataspace::create(res->getTrajectory()), dcpl, lcpl);
      dset_ramo_traj.write(res->getTrajectory());

      auto current_group = part_group.create_group("current");

      auto j = 0;
      for (const auto &cv : res->getCurrent()) {
        std::string ind_sufix =
            std::string(5 - std::to_string(j).length(), '0') +
            std::to_string(j);
        j++;
        auto dset_ramo_current = current_group.create_dataset(
            "current_" + ind_sufix, datatype::create<std::vector<double>>(),
            dataspace::create(cv), dcpl, lcpl);
        dset_ramo_current.write(cv);
      }
    } else {
      usleep(50000);
    }
  }
}

int main(int argc, char *argv[]) {
  using namespace SpaceCharge;
  using namespace hdf5;
  typedef FieldMap<double> fmap;
  typedef std::unique_ptr<fmap> upfmap;
  typedef std::unique_ptr<Track<double>> uptrack;

  Logger::Init();

  auto file = file::open(argv[1], file::AccessFlags::READONLY);
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
  auto nb_part = 1000;
  tbb::concurrent_vector<SpaceCharge::quadv<double>> pos(nb_part);
  for (auto &p : pos) {
    p(0) = 0.0;
    p(1) = pos_x(g1);
    p(2) = pos_y(g1);
    p(3) = pos_z(g1);
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

  std::string filename_tracking = "tracking.h5";
  auto file_tracking =
      file::create(filename_tracking, hdf5::file::AccessFlags::TRUNCATE);

  auto tr_root_group = file_tracking.root();

  auto ind_tr = 0;
  for (auto &track : tracks) {
    std::string ind_tr_sufix =
        std::string(5 - std::to_string(ind_tr).length(), '0') +
        std::to_string(ind_tr);
    ind_tr++;
    node::Group group =
        tr_root_group.create_group("track_" + ind_tr_sufix, lcpl);
    auto dset = group.create_dataset(
        "pos", datatype::create<std::vector<SpaceCharge::quadv<double>>>(),
        dataspace::create(track.getPosVector()), dcpl, lcpl);
    dset.write(track.getPosVector());
    auto dset_speed = group.create_dataset(
        "speed", datatype::create<std::vector<SpaceCharge::quadv<double>>>(),
        dataspace::create(track.getSpeedVector()), dcpl, lcpl);
    dset_speed.write(track.getSpeedVector());
  }

  auto ramofield = readMapFromFile<double>(root_group, "field");
  FieldMapSPS<double> shared = std::move(ramofield);

  std::thread tw_save(save);

  auto ind_part = 0;
  tbb::parallel_for(
      tbb::blocked_range<int>(0, tracks.size()),
      [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i) {
          TrackSPS<double> tr = std::make_shared<Track<double>>(tracks[i]);
          ramoSP ramo = std::make_shared<RamoComputation<double>>(shared, tr);
          ramos.push(std::move(ramo));
        }
      });
  prod = false;
  tw_save.join();
  return 0;
}