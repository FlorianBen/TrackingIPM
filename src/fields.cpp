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

template <class T>
void FieldBunch<T>::addBunch(std::unique_ptr<Bunch<T>> bunch) {
  bunches.push_back(std::move(bunch));
};

template <class T> void FieldBunch<T>::usePeriodicity(bool use) {
  use_periodicity = use;
};

template <class T> T FieldBunch<T>::potentialAt(quadv quad) { return 0.0; };

template <class T> Eigen::Matrix<T, 4, 1> FieldBunch<T>::EfieldAt(quadv quad) {
  Eigen::Matrix<T, 4, 1> E;
  E << 0.0, 0.0, 0.0, 0.0;
  for (auto &bunch : bunches) {
    if (use_periodicity) {
      auto tloc = (quad(0) / cst::sol) + bunch->getBunchPeriod();
      auto nearest =
          std::floor((tloc * 1e9) / (2 * bunch->getBunchPeriod() * 1e9));
      if (nearest > 0) {
        auto shift = 2.0 * nearest * bunch->getBunchPeriod();
        auto shift2 = quad(0) / cst::sol - shift;
        quad(0) = shift2 * cst::sol;
      } else {
      }
      E += bunch->EfieldAt(quad);
      quad(0) = quad(0) - bunch->getBunchPeriod() * cst::sol;
      E += bunch->EfieldAt(quad);
      quad(0) = quad(0) + bunch->getBunchPeriod() * cst::sol;
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