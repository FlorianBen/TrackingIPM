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

namespace SpaceCharge {
template <typename T>
FieldMapSP<T> readMapFromFile(const hdf5::node::Group group, std::string fieldmap_name) {
  auto dataset = group.get_dataset(fieldmap_name);
  hdf5::dataspace::Simple dataspace(dataset.dataspace());

  auto Dimensions = dataspace.current_dimensions();

  auto attr = dataset.attributes;

  quadv<size_t> sizes{0, Dimensions[1], Dimensions[0], Dimensions[2]};
  quadv<T> steps{0., 0., 0, 0.};
  attr["steps"].read(steps);
  quadv<T> offsets{0., 0, 0., 0.};
  attr["offsets"].read(offsets);
  T time = .0;
  attr["time"].read(time);

  FieldMapSP<T> map =
      std::make_unique<FieldMap<T>>(sizes, steps, offsets, time);

  dataset.read(map->getVector());

  return std::move(map);
}

template <typename T> FieldMapSP<T> writeMapFile(const std::string filename) {
  auto file = hdf5::file::open(filename, hdf5::file::AccessFlags::READONLY);
  auto root_group = file.root();
  auto dataset = root_group.get_dataset("field");
  hdf5::dataspace::Simple dataspace(dataset.dataspace());

  auto Dimensions = dataspace.current_dimensions();

  auto attr = dataset.attributes;

  quadv<size_t> sizes{0, Dimensions[1], Dimensions[0], Dimensions[2]};
  quadv<T> steps{0., 0., 0, 0.};
  attr["steps"].read(steps);
  quadv<T> offsets{0., 0, 0., 0.};
  attr["offsets"].read(offsets);
  T time = .0;
  attr["time"].read(time);

  FieldMapSP<T> map =
      std::make_unique<FieldMap<T>>(sizes, steps, offsets, time);

  dataset.read(map->getVector());

  return std::move(map);
}

} // namespace SpaceCharge

#endif