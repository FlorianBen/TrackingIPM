#include "SpaceCharge/field/field_fdm.hpp"

#include "SpaceCharge/core/alogger.hpp"

#include <fstream>
#include <iostream>

namespace SpaceCharge {

FieldFDM::FieldFDM(int nx, int ny, double dx, double dy)
    : nx(nx), ny(ny), size_mat(nx * ny), dx(dx), dy(dy), mat(nx * ny, nx * ny),
      b(nx * ny) {
  Logger::GetLogger()->info("FieldFDM: Matrix ({},{}) created", nx, ny);
}

FieldFDM::~FieldFDM() {}

void FieldFDM::initBoundary() {}

void FieldFDM::initMatrix() {
  std::vector<Tri> coefficients_t;
  auto strip_size = 50;
  for (auto j = 0; j < ny; j++) {
    for (auto i = 0; i < nx; i++) {
      auto indl = gindex(i, j);
      if (j == 0) {
        if ((i < strip_size) || (i > (nx - strip_size))) {
          b(indl) = 0;
          coefficients_t.push_back(Tri(indl, indl, 1));
        } else if (((i > nx / 2 - strip_size) && (i < nx / 2 + strip_size))) {
          b(indl) = 1;
          coefficients_t.push_back(Tri(indl, indl, 1));
        } else if (((i > nx / 4 - strip_size) && (i < nx / 4 + strip_size))) {
          b(indl) = 0;
          coefficients_t.push_back(Tri(indl, indl, 1));
        } else if (((i > 3 * nx / 4 - strip_size) && (i < 3 * nx / 4 + strip_size))) {
          b(indl) = 0;
          coefficients_t.push_back(Tri(indl, indl, 1));
        } else {
          b(indl) = 0;
          coefficients_t.push_back(Tri(indl, gindex(i - 1, j), 1.0 / 2));
          coefficients_t.push_back(Tri(indl, gindex(i + 1, j), 1.0 / 2));
          coefficients_t.push_back(Tri(indl, gindex(i, j + 1), 1));
          coefficients_t.push_back(Tri(indl, indl, -2));
        }
        continue;
      }
      if (j == (ny - 1)) {
        b(indl) = -0;
        coefficients_t.push_back(Tri(indl, indl, 1));
        continue;
      }
      if (i == 0) {
        b(indl) = 0;
        coefficients_t.push_back(Tri(indl, gindex(i, j - 1), 1.0 / 2));
        coefficients_t.push_back(Tri(indl, gindex(i, j + 1), 1.0 / 2));
        coefficients_t.push_back(Tri(indl, gindex(i + 1, j), 1));
        coefficients_t.push_back(Tri(indl, indl, -2));
        continue;
      }
      if (i == (nx - 1)) {
        b(indl) = 0;
        coefficients_t.push_back(Tri(indl, gindex(i, j - 1), 1.0 / 2));
        coefficients_t.push_back(Tri(indl, gindex(i, j + 1), 1.0 / 2));
        coefficients_t.push_back(Tri(indl, gindex(i - 1, j), 1));
        coefficients_t.push_back(Tri(indl, indl, -2));
        continue;
      }
      coefficients_t.push_back(Tri(indl, gindex(i, j - 1), 1));
      coefficients_t.push_back(Tri(indl, gindex(i, j + 1), 1));
      coefficients_t.push_back(Tri(indl, gindex(i - 1, j), 1));
      coefficients_t.push_back(Tri(indl, gindex(i + 1, j), 1));
      coefficients_t.push_back(Tri(indl, indl, -4));
    }
  }

  Logger::GetLogger()->info("Matrix init");
  mat.setFromTriplets(coefficients_t.begin(), coefficients_t.end());

  // Logger::GetLogger()->info("Save fichier");
  // std::ofstream file("mat.txt");
  // if (file.is_open()) {
  //   file << mat << std::endl;
  // }
  // std::cout << b << std::endl;
  // Logger::GetLogger()->info("({}", b);
}

void FieldFDM::solve() {
  Logger::GetLogger()->info("FieldFDM: Solve system on {} thread",
                            Eigen::nbThreads());
  // Eigen::SimplicialCholesky<SpMat> solver(
  //    mat); // performs a Cholesky factorization of A
  // Eigen::VectorXd x = chol.solve(b);
  //, Eigen::RowMajor
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>,
  //                         Eigen::Lower | Eigen::Upper>
  //    solver;
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(1e-6);
  solver.compute(mat);
  x = solver.solve(b);
  Logger::GetLogger()->info("FieldFDM: Solving done in {} iterations, error {}",
                            solver.iterations(), solver.error());
}

void FieldFDM::save(std::string filename) {
  Logger::GetLogger()->info("Save matrix in {}", filename);
  std::ofstream file(filename);
  if (file.is_open()) {
    file << x << std::endl;
  }
}

inline int FieldFDM::gindex(int i, int j) { return i + nx * j; }

} // namespace SpaceCharge
