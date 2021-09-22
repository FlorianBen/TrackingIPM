#ifndef FIELD_MAP_HPP
#define FIELD_MAP_HPP

#include <tbb/concurrent_vector.h>
#include <vector>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/field/fields.hpp"

namespace SpaceCharge {

/**
 * \class FieldMap field_map.hpp
 * \brief FieldMap class describes mapped 3D fields.
 **/
template <typename T> class FieldMap {
protected:
  quadv<size_t> qsize;
  quadv<T> qstep;
  quadv<T> qoffset;
  T time;
  std::vector<quadv<T>> data_;

public:
  /**
   * \brief Default constructor.
   * Construct an empty Track.
   **/
  FieldMap();

  /**
   * \brief Constructor.
   * \param[in] size size vector {x,y,z}.
   * \param[in] step steps vector {x,y,z}.
   * \param[in] offset offset vector {x,y,z}.
   * \param[in] time Time.
   **/
  FieldMap(quadv<size_t> size, quadv<T> step, quadv<T> offset, T time);

  /**
   * \brief Get the total size of the fieldmap.
   * \return Total size.
   **/
  size_t size() const;

  /**
   * \brief Get the dim size of the fieldmap.
   * \return Dim size.
   **/
  quadv<size_t> sizes() const;

  /**
   * \brief Get the steps of the fieldmap.
   * \return Dim size.
   **/
  quadv<T> steps() const;

  /**
   * \brief Get the offsets of the fieldmap.
   * \return Dim size.
   **/
  quadv<T> offsets() const;

  /**
   * \brief Get the x size of the fieldmap.
   * \return x size.
   **/
  size_t nx() const;

  /**
   * \brief Get the y size of the fieldmap.
   * \return y size.
   **/
  size_t ny() const;

  /**
   * \brief Get the z size of the fieldmap.
   * \return z size.
   **/
  size_t nz() const;

  /**
   * \brief Get the raw data vector.
   * \return data pointer.
   **/
  quadv<T> *data();

  /**
   * \brief Get the const raw data vector.
   * \return data pointer.
   **/
  const quadv<T> *data() const;

  /**
   * \brief Operator () read overload.
   * \param[in] x Index x.
   * \param[in] y Index y.
   * \param[in] z Index z.
   */
  quadv<T> operator()(size_t x, size_t y, size_t z) const;

  /**
   * \brief Operator () write overload.
   * \param[in] x Index x.
   * \param[in] y Index y.
   * \param[in] z Index z.
   */
  quadv<T> &operator()(size_t x, size_t y, size_t z);

  void PrintData() const;
}; // namespace SpaceCharge

template <typename T> class FieldMapInterpolate : public FieldMap<T> {

private:
  std::vector<FieldSP<T>> input_fields;



public:
  using FieldMap<T>::FieldMap;

  /**
   * \brief Add a field to the fields manager.
   * \param[in] field Fiel to add.
   */
  void addField(FieldSP<T> &field);

  void computeField();

  void computePot();
};

} // namespace SpaceCharge

#endif