#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <memory>

#include "nanoflann.hpp"
#include "point_cloud.hpp"

namespace SpaceCharge {

/**
 * \class Field fields.hpp
 * \brief Virutal class that represents an EM field.
 **/
template <class T> class EMField {
private:
public:
  /**
   * \brief Constructor.
   * Pure virtual class.
   **/
  EMField();
  virtual ~EMField();

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

  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   * The derived class must implement this function.
   **/
  virtual state_type2<T> EMfieldAt(quadv<T> quad) const = 0;
};

template <typename T> using FieldSP = std::unique_ptr<EMField<T>>;
template <typename T> using FieldSPS = std::shared_ptr<EMField<T>>;

/**
 * \class ConstantField fields.hpp
 * \brief A class that represents a constant EM field.
 **/
template <class T> class ConstantEMField : public EMField<T> {
private:
  state_type2<T> emfield;

public:
  ConstantEMField(state_type2<T> emfield);

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

  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   * The derived class must implement this function.
   **/
  virtual state_type2<T> EMfieldAt(quadv<T> quad) const override;
};

/**
 * \class FieldCOMSOL fields.hpp
 * \brief A class that represents an EM field from a COMSOL file.
 **/
template <class T> class CSVFileEMField : public EMField<T> {
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
  CSVFileEMField();
  virtual ~CSVFileEMField();

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

  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   * The derived class must implement this function.
   **/
  virtual state_type2<T> EMfieldAt(quadv<T> quad) const override;

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
template <class T> class EMFieldsManager : public EMField<T> {
private:
  std::vector<FieldSP<T>> fields;

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

  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   * The derived class must implement this function.
   **/
  virtual state_type2<T> EMfieldAt(quadv<T> quad) const override;
};

}; // namespace SpaceCharge

#endif