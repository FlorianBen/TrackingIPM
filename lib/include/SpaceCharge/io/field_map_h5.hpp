#ifndef FIELD_MAP_H5_HPP
#define FIELD_MAP_H5_HPP

#include <h5cpp/hdf5.hpp>

#include "SpaceCharge/field/field_map.hpp"
#include "SpaceCharge/io/quadv_h5.hpp"

namespace hdf5 {
namespace datatype {

template <typename T> class TypeTrait<SpaceCharge::FieldMap<T>> {
public:
  using Type = SpaceCharge::FieldMap<T>;
  using TypeClass = typename TypeTrait<SpaceCharge::quadv<T>>::TypeClass;

  static TypeClass create(const Type & = Type()) {
    return TypeTrait<SpaceCharge::quadv<T>>::create();
  }
};

} // namespace datatype

namespace dataspace {

template <typename T> class TypeTrait<SpaceCharge::FieldMap<T>> {
public:
  using DataspaceType = Simple;

  static DataspaceType create(const SpaceCharge::FieldMap<T> &value) {
    return Simple(hdf5::Dimensions{value.ny(), value.nx(), value.nz()},
                  hdf5::Dimensions{value.ny(), value.nx(), value.nz()});
  }
  static void *ptr(SpaceCharge::FieldMap<T> &value) {
    return reinterpret_cast<void *>(value.data());
  }

  static const void *cptr(const SpaceCharge::FieldMap<T> &value) {
    return reinterpret_cast<const void *>(value.data());
  }
};

} // namespace dataspace
} // namespace hdf5

#endif