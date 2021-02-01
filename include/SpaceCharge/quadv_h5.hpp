#ifndef QUADV_H5_HPP
#define QUADV_H5_HPP

#include <h5cpp/hdf5.hpp>

namespace hdf5 {
namespace datatype {

template <typename T> class TypeTrait<SpaceCharge::quadv<T>> {
public:
  using Type = SpaceCharge::quadv<T>;
  using TypeClass = Compound;

  static TypeClass create(const Type & = Type()) {
    Compound type = Compound::create(sizeof(Type));
    type.insert("t", 0, TypeTrait<T>::create());
    type.insert("Fx", sizeof(T), TypeTrait<T>::create());
    type.insert("Fy", sizeof(T) * 2, TypeTrait<T>::create());
    type.insert("Fz", sizeof(T) * 3, TypeTrait<T>::create());
    return type;
  }
};

// template <typename T> class TypeTrait<T> {
// public:
//   using Type = T;
//   using TypeClass = Compound;

//   static TypeClass create(const Type & = Type()) {
//     datatype::Compound type =
//         datatype::Compound::create(sizeof(AstroSYCL::core::MonoPixel<T>));
//     type.insert("V", 0, datatype::create<T>());
//     return type;
//   }
// };

} // namespace datatype
} // namespace hdf5

namespace SpaceCharge {} // namespace SpaceCharge

#endif