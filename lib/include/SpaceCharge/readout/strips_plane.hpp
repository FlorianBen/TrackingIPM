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

public:
  /**
   * @brief Construct a new Strips Plane object according to a file.
   *
   * @param filename The text file that contains the strips positions.
   */
  StripsPlane(const std::string filename);

  /**
   * @brief Construct a new Strips Plane object according to a file.
   * This constructor allows to specify the size of the 2D grid using quad
   * vector (t,x,y,z). Only x,y size are revelant.
   * @param filename The text file that contains the strips positions.
   * @param sizes The quadvector size.
   */
  StripsPlane(const std::string filename,
              const SpaceCharge::quadv<size_t> sizes);

  /**
   * @brief Solve the the Ramo problem on the given strip.
   * @param strip Strip to set at 1V.
   */
  void solvePotential(const int strip);

  /**
   * @brief Get the computed potential on the grid.
   * This method will transfert the computed potential into a given FieldMap
   * object.
   * @param fieldmap Input FieldMap object.
   * @param zindex Offset in the FieldMap object.
   */
  void getPotential(FieldMap<double> &fieldmap, int zindex = 0) const;

  /**
   * @brief Get the computed potential on the grid.
   * This method will transfert the computed potential into a given FieldMap
   * object.
   * @param fieldmap Input FieldMap object.
   * @param zindex Offset in the FieldMap object.
   */
  void getField(FieldMap<double> &fieldmap, int zindex = 0) const;

  /**
   * @brief Get the X size of the grid.
   *
   */
  size_t getSizeX() const;
  /**
   * @brief Get the X size of the grid.
   *
   */
  size_t getSizeY() const;

  /**
   * @brief Get the X gap size of the grid.
   *
   */
  double getGapX() const;

  /**
   * @brief Get the X gap size of the grid.
   *
   */
  double getGapY() const;

  /**
   * @brief Get the number of strips.
   * 
   * @return size_t Number of strips
   */
  size_t getNbStrips() const;

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
};

} // namespace SpaceCharge

#endif