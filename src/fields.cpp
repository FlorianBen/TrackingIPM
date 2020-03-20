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
    : use_periodicity(true), local_time(.0){

                             };

template <class T>
void FieldBunch<T>::addBunch(std::unique_ptr<Bunch<T>> bunch) {
  bunches.push_back(std::move(bunch));
};

template <class T> void FieldBunch<T>::usePeriodicity(bool use) {
  use_periodicity = use;
};

template <class T> Eigen::Matrix<T, 4, 1> FieldBunch<T>::EfieldAt(quadv quad) {
  quadv E;
  quadv pos1, pos2, pos3;
  E << 0.0, 0.0, 0.0, 0.0;
  for (auto &bunch : bunches) {
    if (use_periodicity) {
      auto tloc = (quad(0) / cst::sol);
      int rem = (int)std::floor(tloc / (bunch->getBunchPeriod()));
      tloc = tloc - rem * (1 * bunch->getBunchPeriod());
      if (tloc > (bunch->getBunchPeriod() / 2)) {
        quad(0) = (tloc - bunch->getBunchPeriod()) * cst::sol;
      } else {
        quad(0) = tloc * cst::sol;
      }

      E += bunch->EfieldAt(quad);
    } else {
      E += bunch->EfieldAt(quad);
    }
  }

  return E;
}; // namespace SpaceCharge

template <class T>
Eigen::Matrix<T, 4, 1> FieldBunch<T>::MagfieldAt(quadv quad) {

  return quad;
};

template class Field<double>;
template class FieldBunch<double>;
}; // namespace SpaceCharge