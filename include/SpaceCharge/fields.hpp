#ifndef FIELDS_HPP
#define FIELDS_HPP

#include "bunch.hpp"
#include <memory>

namespace SpaceCharge {
template <class T> class Field {
private:
  typedef Eigen::Matrix<T, 4, 1> quadv;

public:
  Field();
  virtual ~Field();

  virtual T potentialAt(quadv quad) = 0;
  virtual quadv EfieldAt(quadv quad) = 0;
  virtual quadv MagfieldAt(quadv quad) = 0;
};

template <class T> class FieldBunch : public Field<T> {
private:
  typedef Eigen::Matrix<T, 4, 1> quadv;
  std::vector<std::unique_ptr<Bunch<T>>> bunches;
  bool use_periodicity;

public:
  FieldBunch();

  void addBunch(std::unique_ptr<Bunch<T>> bunch);
  void usePeriodicity(bool use = true);
  virtual T potentialAt(quadv quad) override;
  virtual quadv EfieldAt(quadv quad) override;
  virtual quadv MagfieldAt(quadv quad) override;
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

  // void addField(Field<T> field);
};

}; // namespace SpaceCharge

#endif