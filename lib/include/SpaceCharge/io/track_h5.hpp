#ifndef TRACK_H5_HPP
#define TRACK_H5_HPP

#include <h5cpp/hdf5.hpp>

#include "SpaceCharge/io/quadv_h5.hpp"
#include "SpaceCharge/track/track.hpp"

namespace hdf5 {

namespace dataspace {

template <typename T> class TypeTrait<std::vector<SpaceCharge::quadv<T>>> {
public:
  using DataspaceType = Simple;

  static DataspaceType create(const std::vector<SpaceCharge::quadv<T>> &value) {
    return Simple(hdf5::Dimensions{value.size()});
  }
  static void *ptr(std::vector<SpaceCharge::quadv<T>> &value) {
    return reinterpret_cast<void *>(value.data());
  }

  static const void *cptr(const std::vector<SpaceCharge::quadv<T>> &value) {
    return reinterpret_cast<const void *>(value.data());
  }
};

} // namespace dataspace
} // namespace hdf5

#endif