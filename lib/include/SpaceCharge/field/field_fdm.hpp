#ifndef FIELD_FDM_HPP
#define FIELD_FDM_HPP

#include <Eigen/Eigen>

namespace SpaceCharge {
class FieldFDM {
  typedef Eigen::SparseMatrix<double>
      SpMat; // declares a column-major sparse matrix type of double
  typedef Eigen::Triplet<double> Tri;

private:
  /* data */
  int nx;
  int ny;
  double dx;
  double dy;
  double t;

  SpMat mat;


public:
  FieldFDM(int nx, int ny, double dx, double dy);
  ~FieldFDM();

  void init_boundary();
  void init_rhs();
  void init_matrix();
  void solve();

};
} // namespace SpaceCharge

#endif