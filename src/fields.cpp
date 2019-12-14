#include "SpaceCharge/fields.hpp"

#include <iostream>

namespace SpaceCharge {
template <class T>
Field<T>::Field(){

};

template <class T>
Field<T>::~Field(){

};

template <class T>
FieldBunch<T>::FieldBunch()
    : use_periodicity(true){

      };

//template <class T> void FieldBunch<T>::addBunch(Bunch<T> bunch) {
  //bunches.push_back(bunch);
//};

template <class T> void FieldBunch<T>::usePeriodicity(bool use) {
  use_periodicity = use;
};

template <class T> T FieldBunch<T>::potentialAt(quadv quad) {
  return 0.0;
};

template <class T> Eigen::Matrix<T, 4, 1> FieldBunch<T>::EfieldAt(quadv quad) {
  return quad;
};

template <class T>
Eigen::Matrix<T, 4, 1> FieldBunch<T>::MagfieldAt(quadv quad) {

  return quad;
};

template class Field<double>;
template class FieldBunch<double>;

}; // namespace SpaceCharge