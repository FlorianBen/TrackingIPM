#pragma once

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
  /** Quadrivector time - position*/
  typedef Eigen::Matrix<T, 4, 1> quadv;

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
  virtual quadv EfieldAt(quadv quad) = 0;
  /**
   * \brief Calculate the Magnetic field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   * The derived class must implement this function.
   **/
  virtual quadv MagfieldAt(quadv quad) = 0;
};

/**
 * \class ConstantField fields.hpp
 * \brief A class that represents a constant EM field.
 **/
template <class T> class ConstantField : public Field<T> {
private:
  /** Quadrivector time - position*/
  typedef Eigen::Matrix<T, 4, 1> quadv;
  quadv field;

public:
  ConstantField(quadv field);
  virtual ~ConstantField();

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   * The returned value is constant.
   **/
  virtual quadv EfieldAt(quadv quad) override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   * The returned value is constant.
   **/
  virtual quadv MagfieldAt(quadv quad) override;
};

/**
 * \class FieldBunch fields.hpp
 * \brief A class that represents an EM field created by particle bunches.
 * Several bunch can be added. The resulting EM field is the sum of the
 * contribution of each bunch.
 **/
template <class T> class FieldBunch : public Field<T> {
private:
  /** Quadrivector time - position*/
  typedef Eigen::Matrix<T, 4, 1> quadv;
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
  virtual quadv EfieldAt(quadv quad) override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   **/
  virtual quadv MagfieldAt(quadv quad) override;
};

/**
 * \class FieldCOMSOL fields.hpp
 * \brief A class that represents an EM field from a COMSOL file.
 **/
template <class T> class FieldCOMSOL : public Field<T> {
private:
  /** Quadrivector time - position*/
  typedef Eigen::Matrix<T, 4, 1> quadv;
  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<T, PointCloud<T>>, PointCloud<T>, 4 /* dim */
      >
      kd_tree_nanoflann;

  PointCloud<T> pointcloud_posE;
  PointCloud<T> fieldE;
  kd_tree_nanoflann *index;

public:
  FieldCOMSOL();
  virtual ~FieldCOMSOL();

  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   **/
  virtual quadv EfieldAt(quadv quad) override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   **/
  virtual quadv MagfieldAt(quadv quad) override;

  size_t interpolateNN(const quadv pos, quadv &fieldv, size_t size = 7) const;
  size_t interpolateNN2(const quadv pos, quadv &fieldv, size_t size = 7) const;

  void interpolateRBF(const quadv pos, quadv &fieldv, const int nbNN,
                      const float order,
                      const std::function<T(const T, const T)> kernel) const;

  void loadEfield(const std::string filename, const quadv offset,
                  const double scale = 1.0);
  void loadBfield(const std::string filename, const quadv offset,
                  const double scale = 1.0);

  void create_Eindex(const int leaf_size = 10);
};

/**
 * \class Fields fields.hpp
 * \brief A class that represents an EM field created by several fields.
 * The resulting EM field is the sum of the
 * contribution of each field.
 **/
template <class T> class Fields {
private:
  std::vector<Field<T>> fields;

public:
  Fields();

  // void addField(Field<T> field);
};

}; // namespace SpaceCharge