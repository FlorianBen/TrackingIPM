#ifndef FIELD_MAP_HPP
#define FIELD_MAP_HPP

#include <tbb/concurrent_vector.h>
#include <vector>

#include "SpaceCharge/alogger.hpp"
#include "SpaceCharge/fields.hpp"

namespace SpaceCharge {

template <typename T> class FieldMap {
private:
  quadv<size_t> qsize;
  quadv<T> qstep;
  quadv<T> qoffset;

  tbb::concurrent_vector<quadv<T>> data_;
  std::vector<quadv<T>> data2_;
  tbb::concurrent_vector<T> data_pot;
  std::vector<FieldSP<T>> input_fields;

public:
  FieldMap();
  FieldMap(quadv<size_t> size, quadv<T> step, quadv<T> offset);

  const quadv<T> &operator()(size_t i, size_t j, size_t k) const;

  quadv<T> &operator()(size_t i, size_t j, size_t k);

  size_t size() const;

  size_t nx() const;

  size_t ny() const;

  size_t nz() const;

  quadv<T> *data();

  const quadv<T> *data() const;

  /**
   * \brief Add a field to the fields manager.
   * \param[in] field Fiel to add.
   */
  void addField(FieldSP<T> &field);

  void computeField();

  void computePot();

  void PrintData() const;
}; // namespace SpaceCharge

} // namespace SpaceCharge

#endif