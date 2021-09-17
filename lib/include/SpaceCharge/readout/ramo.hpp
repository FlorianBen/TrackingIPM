#ifndef RAMO_HPP
#define RAMO_HPP

#include <memory>
#include <vector>

#include "SpaceCharge/core/nanoflann.hpp"
#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/track/track.hpp"

namespace SpaceCharge {

/**
 * @brief This class computes the signal induced on conductive strips according
 * to the Ramo Shockley.
 *
 * @tparam T
 */
template <typename T> class RamoComputation {
  typedef std::unique_ptr<FieldMap<T>> fmap;
  typedef std::unique_ptr<Track<T>> ftrack;
  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<T, PointCloud<T>>, PointCloud<T>, 4 /* dim */
      >
      kd_tree_nanoflann;
  typedef std::vector<std::vector<T>> testssfd;

private:
  int size;
  fmap map;
  ftrack track;
  PointCloud<T> pointcloud_pos;
  kd_tree_nanoflann *index;
  state_type2<int> coordinates;
  state_type2<int> res_coordinates;

  state_type2<T> tracjectory;
  std::vector<std::vector<T>> current;

  void createPosCloud();

  void createIndex(const int leaf_size = 20);

  void computeTrajectory();

  void computeCurrent();

  void computeStrips();

public:
  RamoComputation(fmap &map, ftrack &track);

  ~RamoComputation();

  const state_type2<int>& getCoordinates() const;

  const state_type2<T>& getTrajectory() const;

  const std::vector<std::vector<T>>& getCurrent() const;
};
} // namespace SpaceCharge
#endif