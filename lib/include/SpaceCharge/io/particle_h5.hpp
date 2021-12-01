#ifndef PARTICLE_H5_HPP
#define PARTICLE_H5_HPP

#include <h5cpp/hdf5.hpp>

// namespace hdf5 {
// namespace datatype {

// template <typename T> class TypeTrait<SpaceCharge::<T>> {
// public:
//   using Type = SpaceCharge::quadv<T>;
//   using TypeClass = Compound;

//   static TypeClass create(const Type & = Type()) {
//     Compound type = Compound::create(sizeof(Type));
//     type.insert("t", 0, TypeTrait<T>::create());
//     type.insert("x", sizeof(T), TypeTrait<T>::create());
//     type.insert("y", sizeof(T) * 2, TypeTrait<T>::create());
//     type.insert("z", sizeof(T) * 3, TypeTrait<T>::create());
//     return type;
//   }
// };

} // namespace datatype
} // namespace hdf5
