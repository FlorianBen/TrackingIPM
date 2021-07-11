#ifndef FIELD_FDM_HPP
#define FIELD_FDM_HPP

#include <Eigen/Eigen>

namespace SpaceCharge {
class FieldFDM {
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor>//, Eigen::RowMajor>
      SpMat; // declares a column-major sparse matrix type of double
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

public:
  FieldFDM(int nx, int ny, double dx, double dy);
  ~FieldFDM();

  void initBoundary();
  void init_rhs();
  void initMatrix();
  void solve();

  void save(std::string filename);

  int gindex(int i, int j);
};
} // namespace SpaceCharge

#endif