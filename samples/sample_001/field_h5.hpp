#pragma once

#include <h5cpp/hdf5.hpp>
#include <tbb/concurrent_vector.h>

#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/field/fields.hpp"
#include "SpaceCharge/core/particle.hpp"

#include "vector.hpp"

namespace hdf5 {
namespace datatype {

//!
//! \brief datatype trait for a 3D vector
//!
template <typename T> class TypeTrait<Vector<T>> {
public:
  using Type = Vector<T>;
  using TypeClass = Compound;

  static TypeClass create(const Type & = Type()) {
    Compound type = Compound::create(sizeof(Vector<T>));
    type.insert("x", 0, TypeTrait<T>::create());
    type.insert("y", sizeof(T), TypeTrait<T>::create());
    type.insert("z", sizeof(T) * 2, TypeTrait<T>::create());
    return type;
  }
};

} // namespace datatype
} // namespace hdf5

using FloatVector = Vector<float>;
using Field = tbb::concurrent_vector<FloatVector>;

class VectorAppender {
private:
  hdf5::node::Dataset dataset_;
  hdf5::dataspace::Hyperslab selection_;
  std::string log_prefix_;

public:
  VectorAppender() = delete;
  VectorAppender(hdf5::node::Dataset dataset, const std::string &log_prefix)
      : dataset_(dataset),
        selection_(hdf5::dataspace::Hyperslab({0}, {1}, {1}, {1})),
        log_prefix_(log_prefix) {
    hdf5::dataspace::Simple dataspace = dataset_.dataspace();
    selection_.offset(0, dataspace.current_dimensions()[0]);
  }

  void operator()(const FloatVector vector) {
    dataset_.extent(0, 1);                            // extend the dataset
    dataset_.write(vector, selection_);               // write the data
    selection_.offset(0, selection_.offset()[0] + 1); // extend the selection
  }
};

// hdf5::node::Dataset
hdf5::node::Dataset
create_vector_dataset(const std::string &name,
                      const hdf5::node::Group &parent_group) {
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  dcpl.chunk(hdf5::Dimensions{256});

  auto datatype = hdf5::datatype::TypeTrait<FloatVector>::create();
  hdf5::dataspace::Simple dataspace({0}, {hdf5::dataspace::Simple::UNLIMITED});
  return hdf5::node::Dataset(parent_group, name, datatype, dataspace, lcpl,
                             dcpl);
}