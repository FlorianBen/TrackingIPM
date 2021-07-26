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

void FieldFDM::readStrips(const std::string filepath) {
  Logger::GetLogger()->info("FieldFDM: Load strip distribution from file.");
  std::ifstream infile(filepath);
  double a, b;
  while (infile >> a >> b) {
    std::pair<double, double> pair_corr{a * 1e-3 + 0.04, b * 1e-3 + 0.04};
    strips_pairs.push_back(pair_corr);
  }
  for (auto pair : strips_pairs) {
    int r1 = std::round(pair.first / dx);
    int r2 = std::round(pair.second / dx);
    // Logger::GetLogger()->info("{} <> {}", r1, r2);
  }
}

bool FieldFDM::isStrips(const int i, int &stripnb) {
  auto ind = 0;
  for (auto pair : strips_pairs) {
    int r1 = std::round(pair.first / dx);
    int r2 = std::round(pair.second / dx);
    if (r1 <= i && i <= r2) {
      stripnb = ind;
      return true;
    }
    ind++;
  }
  return false;
}

void FieldFDM::initMatrix() {
  std::vector<Tri> coefficients_t;
  auto strip_size = 50;
  auto strip_on = 9;
  auto strip_nb = 0;
  for (auto j = 0; j < ny; j++) {
    for (auto i = 0; i < nx; i++) {
      auto indl = gindex(i, j);
      if (j == 0) {
        if (isStrips(i, strip_nb)) {
          if (strip_on == strip_nb) {
            b(indl) = 1;

          } else {
            b(indl) = 0;
          }

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

  Logger::GetLogger()->info("FieldFDM: Matrix initialization");
  mat.setFromTriplets(coefficients_t.begin(), coefficients_t.end());
}

void FieldFDM::solve() {
  Logger::GetLogger()->info("FieldFDM: Solve system on {} thread",
                            Eigen::nbThreads());
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(1e-9);
  solver.compute(mat);
  x = solver.solve(b);
  Logger::GetLogger()->info("FieldFDM: Solving done in {} iterations, error {}",
                            solver.iterations(), solver.error());
}

void FieldFDM::save(std::string filename) {
  Logger::GetLogger()->info("FieldFDM: Save matrix in {}", filename);
  std::ofstream file(filename);
  if (file.is_open()) {
    file << x << std::endl;
  }
}

inline int FieldFDM::gindex(int i, int j) { return i + nx * j; }

} // namespace SpaceCharge
