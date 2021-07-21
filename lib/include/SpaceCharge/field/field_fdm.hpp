#ifndef FIELD_FDM_HPP
#define FIELD_FDM_HPP

#include <Eigen/Eigen>

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

  SpMat mat;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  int gindex(int i, int j);

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

  /**
   * @brief Export field to file.
   * 
   * @param filename Filename.
   */
  void save(std::string filename);
};
} // namespace SpaceCharge

#endif