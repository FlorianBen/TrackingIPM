#include "SpaceCharge/fields.hpp"

namespace SpaceCharge {
template <class T>
Field<T>::Field(){

};

template class Field<double>;
template class Field<float>;
}; // namespace SpaceCharge