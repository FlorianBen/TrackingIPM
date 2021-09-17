#include "SpaceCharge/readout/ramo.hpp"
#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/nanoflann.hpp"
#include "SpaceCharge/core/point_cloud.hpp"

namespace SpaceCharge {

// TODO: Add more generic version.

template <typename T>
RamoComputation<T>::RamoComputation(fmap &map, ftrack &track)
    : size(100), map(std::move(map)), track(std::move(track)) {
  SC_INFO("RamoComputation:: Constructor");
  computeTrajectory();
}

template <typename T> void RamoComputation<T>::computeTrajectory() {

  auto pos = track->getPosVector().back();
  auto vit = track->getSpeedVector().back();

  auto d_max = map->sizes()(2) * map->steps()(2);
  auto t_max = d_max / vit(1);

  for (int y = 0; y < map->sizes()(2); y++) {
    for (int x = 0; x < map->sizes()(1); x++) {
      quadv<T> tmp_pos(0.0, x * map->steps()(1), y * map->steps()(2), 0.0);
    }
  }

  for (auto i = 0; i < size; i++) {
    quadv<T> tmp_pos{0.0, (1.0 * i / size) * t_max * vit(1),
                     (1.0 * i / size) * t_max * vit(2), 0.0};
    tracjectory.push_back(tmp_pos);
  }
}

template <typename T> RamoComputation<T>::~RamoComputation() {}

template class RamoComputation<double>;

} // namespace SpaceCharge