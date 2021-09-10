#ifndef RAMO_HPP
#define RAMO_HPP

#include <memory>

#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/track/track.hpp"

namespace SpaceCharge {

template <typename T> class RamoComputation {
  typedef std::unique_ptr<FieldMap<T>> fmap;
  typedef std::unique_ptr<Track<T>> ftrack;

public:
    RamoComputation();

private:
  fmap map;
  ftrack track;


  void computeCurrent();

  void computeStrips();
};
} // namespace SpaceCharge
#endif