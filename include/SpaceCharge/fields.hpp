#ifndef FIELDS_HPP
#define FIELDS_HPP

#include "bunch.hpp"

namespace SpaceCharge {
template <class T> class Field {
private:

public:
  Field();
};

template <class T> class FieldBunch {
private:
  std::vector<Bunch<T>> bunches;

public:
  FieldBunch();
};

template <class T> class FieldCOMSOL {
private:

public:
  FieldCOMSOL();
};

template <class T> class Fields {
private:
  std::vector<Field<T>> fields;

public:
  Fields();

  void addField();
};

}; // namespace SpaceCharge

#endif