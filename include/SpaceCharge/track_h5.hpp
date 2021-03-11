#ifndef TRACK_H5_HPP
#define TRACK_H5_HPP

#include <h5cpp/hdf5.hpp>

#include "SpaceCharge/quadv_h5.hpp"
#include "SpaceCharge/track.hpp"

namespace hdf5 {
namespace datatype {

template <typename T> class TypeTrait<SpaceCharge::Track<T>> {
public:
  using Type = SpaceCharge::Track<T>;
  using TypeClass = typename TypeTrait<SpaceCharge::quadv<T>>::TypeClass;

  static TypeClass create(const Type & = Type()) {
    return TypeTrait<SpaceCharge::quadv<T>>::create();
  }
};

} // namespace datatype

namespace dataspace {

template <typename T> class TypeTrait<SpaceCharge::Track<T>> {
public:
  using DataspaceType = Simple;

  static DataspaceType create(const SpaceCharge::Track<T> &value) {
    return Simple(hdf5::Dimensions{value.size()});
  }
  static void *ptr(SpaceCharge::Track<T> &value) {
    return reinterpret_cast<void *>(value.data());
  }

  static const void *cptr(const SpaceCharge::Track<T> &value) {
    return reinterpret_cast<const void *>(value.data());
  }
};

} // namespace dataspace
} // namespace hdf5

// namespace SpaceCharge {
// template <typename T> class VectorAppender {
// private:
//   hdf5::node::Dataset dataset_;
//   hdf5::dataspace::Hyperslab selection_;
//   std::string log_prefix_;

// public:
//   VectorAppender() = delete;
//   VectorAppender(hdf5::node::Dataset dataset, const std::string &log_prefix)
//       : dataset_(dataset),
//         selection_(hdf5::dataspace::Hyperslab({0}, {1}, {1}, {1})),
//         log_prefix_(log_prefix) {
//     hdf5::dataspace::Simple dataspace = dataset_.dataspace();
//     selection_.offset(0, dataspace.current_dimensions()[0]);
//   }

//   void operator()(const SpaceCharge::quadv<T> vector) {
//     dataset_.extent(0, 1);                            // extend the dataset
//     dataset_.write(vector, selection_);               // write the data
//     selection_.offset(0, selection_.offset()[0] + 1); // extend the selection
//   }
// };

// // // hdf5::node::Dataset
// template <typename T>
// hdf5::node::Dataset
// create_vector_dataset(const std::string &name,
//                       const hdf5::node::Group &parent_group) {
//   hdf5::property::LinkCreationList lcpl;
//   hdf5::property::DatasetCreationList dcpl;
//   dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
//   dcpl.chunk(hdf5::Dimensions{256});

//   auto datatype = hdf5::datatype::TypeTrait<SpaceCharge::quadv<T>>::create();
//   hdf5::dataspace::Simple dataspace({0},
//   {hdf5::dataspace::Simple::UNLIMITED}); return
//   hdf5::node::Dataset(parent_group, name, datatype, dataspace, lcpl,
//                              dcpl);
// }
// } // namespace SpaceCharge

#endif