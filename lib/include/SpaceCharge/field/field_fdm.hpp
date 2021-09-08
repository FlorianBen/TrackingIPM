#ifndef FIELD_FDM_HPP
#define FIELD_FDM_HPP

#include <Eigen/Eigen>

#include "SpaceCharge/field/field_map.hpp"

namespace SpaceCharge {
/**
 * @brief
 *
 */
class FieldFDM {
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
  typedef Eigen::Triplet<double> Tri;

private:
  /* data */
  int nx;
  int ny;
  int size_mat;
  double dx;
  double dy;
  double t;

  int total_strip;

  SpMat mat;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  std::vector<std::pair<double, double>> strips_pairs;

  int gindex(int i, int j) const;

  void initStrips();

  bool isStrips(const int i, int &stripnb);

public:
  /**
   * @brief Construct a new Field FDM object.
   *
   * @param nx Size of the x axis.
   * @param ny Size of the y axis.
   * @param dx Delta between two x points.
   * @param dy Delta between two y points.
   */
  FieldFDM(int nx, int ny, double dx, double dy);
  ~FieldFDM();

  void readStrips(const std::string filepath);

  /**
   * @brief Initialize the FDM matrix.
   *
   */
  void initMatrix();
  /**
   * @brief Solve the the FDM problem.
   *
   */
  void solve();

  void fillFieldMap(FieldMap<double> &fieldmap) const;

  void fillFieldMap2(FieldMap<double> &fieldmap) const;

  /**
   * @brief Export field to file.
   *
   * @param filename Filename.
   */
  void save(std::string filename);
};
} // namespace SpaceCharge

#endif