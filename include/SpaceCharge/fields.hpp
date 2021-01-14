#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <memory>

#include "bunch.hpp"
#include "nanoflann.hpp"
#include "point_cloud.hpp"

namespace SpaceCharge {

/**
 * \class Field fields.hpp
 * \brief Virutal class that represents an EM field.
 **/
template <class T> class Field {
private:
public:
  /**
   * \brief Constructor.
   * Pure virtual class.
   **/
  Field();
  virtual ~Field();

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   * The derived class must implement this function.
   **/
  virtual quadv<T> EfieldAt(quadv<T> quad) const = 0;
  /**
   * \brief Calculate the Magnetic field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   * The derived class must implement this function.
   **/
  virtual quadv<T> MagfieldAt(quadv<T> quad) const = 0;
};

template <typename T> using FieldSP = std::unique_ptr<Field<T>>;
template <typename T> using FieldSPS = std::shared_ptr<Field<T>>;

/**
 * \class ConstantField fields.hpp
 * \brief A class that represents a constant EM field.
 **/
template <class T> class ConstantField : public Field<T> {
private:
  quadv<T> field;

public:
  ConstantField(quadv<T> field);

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   * The returned value is constant.
   **/
  virtual quadv<T> EfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   * The returned value is constant.
   **/
  virtual quadv<T> MagfieldAt(quadv<T> quad) const override;
};

/**
 * \class FieldBunch fields.hpp
 * \brief A class that represents an EM field created by particle bunches.
 * Several bunch can be added. The resulting EM field is the sum of the
 * contribution of each bunch.
 **/
template <class T> class FieldBunch : public Field<T> {
private:
  std::vector<std::unique_ptr<Bunch<T>>> bunches;
  bool use_periodicity;
  T local_time;

public:
  FieldBunch();

  /**
   * \brief Add a bunch as source of EM field.
   * \param[in] bunch Pointer to a bunch.
   **/
  void addBunch(std::unique_ptr<Bunch<T>> bunch);
  /**
   * \brief Use the period of the bunch.
   * \param[in] use Use periodicity if True.
   **/
  void usePeriodicity(bool use = true);
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   **/
  virtual quadv<T> EfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   **/
  virtual quadv<T> MagfieldAt(quadv<T> quad) const override;
};

/**
 * \class FieldCOMSOL fields.hpp
 * \brief A class that represents an EM field from a COMSOL file.
 **/
template <class T> class FieldCOMSOL : public Field<T> {
private:
  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<T, PointCloud<T>>, PointCloud<T>, 4 /* dim */
      >
      kd_tree_nanoflann;

  PointCloud<T> pointcloud_posE;
  PointCloud<T> fieldE;
  kd_tree_nanoflann *index;

  static T kernel_exp(T distance, T n) { return exp(-pow(n * distance, 2)); }

public:
  FieldCOMSOL();
  virtual ~FieldCOMSOL();

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   **/
  virtual quadv<T> EfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   **/
  virtual quadv<T> MagfieldAt(quadv<T> quad) const override;

  size_t interpolateNN(const quadv<T> pos, quadv<T> &fieldv,
                       size_t size = 7) const;
  size_t interpolateNN2(const quadv<T> pos, quadv<T> &fieldv,
                        size_t size = 7) const;

  void interpolateRBF(const quadv<T> pos, quadv<T> &fieldv, const int nbNN,
                      const float order,
                      const std::function<T(const T, const T)> kernel) const;

  void loadEfield(const std::string filename, const quadv<T> offset,
                  const double scale = 1.0, const int leaf_size = 20);
  void loadBfield(const std::string filename, const quadv<T> offset,
                  const double scale = 1.0, const int leaf_size = 20);

  void create_Eindex(const int leaf_size = 20);
  void create_Bindex(const int leaf_size = 20);
};

/**
 * \class Fields fields.hpp
 * \brief A class that represents an EM field created by several fields.
 * The resulting EM field is the sum of the
 * contribution of each field.
 **/
template <class T> class EMFieldsManager : public Field {
private:
  std::vector<FieldSP<T>> Efields;

public:
  EMFieldsManager();

  /**
   * \brief Add a field to the fields manager.
   * \param[in] field Fiel to add.
   */
  void addField(FieldSP<T> &field);

  /**
   * \brief Remove the field n.
   * \param[in] index Index.
   */
  void removeField(std::size_t index);

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   **/
  virtual quadv<T> EfieldAt(quadv<T> quad) const override;

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   **/
  virtual quadv<T> MagfieldAt(quadv<T> quad) const override;
};

}; // namespace SpaceCharge

#endif