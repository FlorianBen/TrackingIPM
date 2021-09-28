#include "SpaceCharge/readout/ramo.hpp"
#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/definitions.hpp"
#include "SpaceCharge/core/nanoflann.hpp"
#include "SpaceCharge/core/point_cloud.hpp"

namespace SpaceCharge {

// TODO: Add more generic version.
// TODO: Select the particle direction.

template <typename T>
RamoComputation<T>::RamoComputation(fmap &map, ftrack &track)
    : size(100), map(map), track(track) {
  createPosCloud();
  computeTrajectory();
  computeCurrent();
}

template <typename T> RamoComputation<T>::~RamoComputation() { delete index; }

template <typename T>
const state_type2<T> &RamoComputation<T>::getTrajectory() const {
  return tracjectory;
}

template <typename T>
const state_type2<int> &RamoComputation<T>::getCoordinates() const {
  return res_coordinates;
}

template <typename T>
const std::vector<std::vector<T>> &RamoComputation<T>::getCurrent() const {
  return current;
}

template <typename T> void RamoComputation<T>::computeTrajectory() {

  auto pos = track->getPosVector().back();

  vit = track->getSpeedVector().back();

  auto d_max = map->sizes()(2) * map->steps()(2);
  auto t_max = d_max / vit(2);

  for (auto i = 0; i < size; i++) {
    quadv<T> tmp_pos{0.0, (1.0 * i / size) * t_max * vit(1) + pos(1),
                     (1.0 * i / size) * t_max * vit(2), 0.0};
    tracjectory.push_back(tmp_pos);
  }
}

template <typename T> void RamoComputation<T>::createPosCloud() {
  for (int y = 0; y < map->sizes()(2); y++) {
    for (int x = 0; x < map->sizes()(1); x++) {
      quadv<T> tmp_pos(0.0, x * map->steps()(1), y * map->steps()(2), 0.0);
      pointcloud_pos.pts.push_back(tmp_pos - map->offsets());
      coordinates.push_back(quadv<int>(0, x, y, 0));
    }
  }
  createIndex();
}

template <typename T>
void RamoComputation<T>::createIndex(const int leaf_size) {
  index = new kd_tree_nanoflann(
      3, pointcloud_pos, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size));
  index->buildIndex();
}

template <typename T> void RamoComputation<T>::computeCurrent() {
  quadv<T> fieldv;
  for (auto i = 0; i < map->sizes()(3); i++) {
    std::vector<T> temp_current;
    for (const auto &pos : tracjectory) {
      T query_pt[4] = {pos[0], pos[1], pos[2], pos[3]};
      size_t nb_result = size;
      std::vector<size_t> ret_index(nb_result);
      std::vector<T> out_dist_sqr(nb_result);
      nb_result = index->knnSearch(&query_pt[0], nb_result, &ret_index[0],
                                   &out_dist_sqr[0]);
      res_coordinates.push_back(coordinates[ret_index[0]]);
      auto coor = coordinates[ret_index[0]];
      fieldv = map->operator()(coor(1), coor(2), i);
      auto temp = cst::e * scalar_prod(fieldv, vit);
      temp_current.push_back(temp);
    }
    current.push_back(temp_current);
  }
}

template class RamoComputation<double>;

} // namespace SpaceCharge