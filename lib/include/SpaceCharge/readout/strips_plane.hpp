#ifndef STRIPS_PLANE_HPP
#define STRIPS_PLANE_HPP

#include <Eigen/Eigen>
#include <string>
#include <vector>

#include "SpaceCharge/field/field_map.hpp"

namespace SpaceCharge {

/**
 * @brief This class describes a 2D strips PCB plane.
 * This class encapsulates the data for the computation of the sensitivity
 * fields. It requires an input text file that contains the pair of positions of
 * each strips in ascending order. The first and last pair are considered as
 * fixed ground not strip.
 */
class StripsPlane {
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
  typedef Eigen::Triplet<double> Tri;

private:
  int nx;
  int ny;
  int size_mat;

  double dx;
  double dy;

  double gap;

  double min_dx;
  double pcb_size;
  int total_strips;

  std::vector<std::pair<double, double>> strips_pairs;

  SpMat mat;
  Eigen::VectorXd b;
  Eigen::VectorXd x;

  void readFile(const std::string filename);

  void computeBestSize();

  void initMatrix(const int strip);
  void solve();
  bool isStrips(const int i, int &stripnb) const;
  int gindex(int i, int j) const;

public:
  /**
   * @brief Construct a new Strips Plane object according to a file.
   *
   * @param filename The text file that contains the strips positions.
   */
  StripsPlane(const std::string filename);

  /**
   * @brief Solve the the FDM problem.
   *
   */
  void solvePotential(const int strip);

  void getPotential(FieldMap<double> &fieldmap) const;

  void getField(FieldMap<double> &fieldmap) const;
};

} // namespace SpaceCharge

#endif