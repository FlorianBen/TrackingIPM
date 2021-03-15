#ifndef FIELD_MAP_HPP
#define FIELD_MAP_HPP

#include <tbb/concurrent_vector.h>
#include <vector>

#include "SpaceCharge/alogger.hpp"
#include "SpaceCharge/fields.hpp"

namespace SpaceCharge {

/**
 * \class FieldMap field_map.hpp
 * \brief FieldMap class describes mapped 3D fields.
 **/
template <typename T> class FieldMap {
private:
  quadv<size_t> qsize;
  quadv<T> qstep;
  quadv<T> qoffset;
  T time;

  std::vector<quadv<T>> data_;
  std::vector<FieldSP<T>> input_fields;

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