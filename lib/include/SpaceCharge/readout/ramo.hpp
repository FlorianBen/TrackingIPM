#ifndef RAMO_HPP
#define RAMO_HPP

#include <memory>

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

public:
  RamoComputation(fmap& map, ftrack& track);

  ~RamoComputation();

private:
  int size;
  fmap map;
  ftrack track;
  PointCloud<T> pointcloud_pos;
  state_type2<T> tracjectory;

  void computeTrajectory();

  void computeCurrent();

  void computeStrips();
};
} // namespace SpaceCharge
#endif