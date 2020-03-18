#pragma once

#include <array>
#include <iostream>

//!
//! \brief custom vector type
//!
template<typename T>
struct Vector : public std::array<T,3>
{
  using std::array<T,3>::array;
  using std::array<T,3>::operator=;

};

template<typename T>
std::ostream &operator<<(std::ostream &stream,const Vector<T> vector)
{
  stream<<"("<<vector[0]<<","<<vector[1]<<","<<vector[2]<<")";
  return stream;
}