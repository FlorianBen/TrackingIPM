
#include <fstream>
#include <iostream>
#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/field/field_fdm.hpp"

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
  total_strip = strips_pairs.size();
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
  auto strip_on = 1;
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
  solver.setTolerance(1e-7);
  solver.compute(mat);
  x = solver.solve(b);
  Logger::GetLogger()->info("FieldFDM: Solving done in {} iterations, error {}",
                            solver.iterations(), solver.error());
}

void FieldFDM::fillFieldMap(FieldMap<double> &fieldmap) const {
  Logger::GetLogger()->info("FieldFDM: Fill a field map");
  tbb::parallel_for(
      tbb::blocked_range3d<int>(0, 1, 0, ny, 0, nx),
      [&](const tbb::blocked_range3d<int> &r) {
        for (int y = r.rows().begin(), y_end = r.rows().end(); y < y_end; y++) {
          for (int x = r.cols().begin(), x_end = r.cols().end(); x < x_end;
               x++) {
            for (int z = r.pages().begin(), z_end = r.pages().end(); z < z_end;
                 z++) {
              quadv<double> test{0.0, this->x(x + nx * y), 0.0, 0.0};
              fieldmap(x, y, z) = test;
            }
          }
        }
      });
}

void FieldFDM::fillFieldMap2(FieldMap<double> &fieldmap) const {
  Logger::GetLogger()->info("FieldFDM: Fill a field map");
  tbb::parallel_for(
      tbb::blocked_range3d<int>(0, 1, 0, ny, 0, nx),
      [&](const tbb::blocked_range3d<int> &r) {
        for (int j = r.rows().begin(), y_end = r.rows().end(); j < y_end; j++) {
          for (int i = r.cols().begin(), x_end = r.cols().end(); i < x_end;
               i++) {
            for (int z = r.pages().begin(), z_end = r.pages().end(); z < z_end;
                 z++) {
              auto gid = gindex(i, j);
              auto Ex = 0.0;
              auto Ey = 0.0;
              auto Ez = 0.0;
              if (i == 0) {
                Ex = (this->x(gindex(i + 1, j)) - this->x(gid)) / (dx);
              } else if (i == (x_end - 1)) {
                Ex = (this->x(gid) - this->x(gindex(i - 1, j))) / (dx);
              } else {
                Ex = (this->x(gindex(i + 1, j)) - this->x(gindex(i - 1, j))) /
                     (2 * dx);
              }
              if (j == 0) {
                Ey = (this->x(gindex(i, j + 1)) - this->x(gid)) / dx;
              } else if (j == (y_end - 1)) {
                Ey = (this->x(gid) - this->x(gindex(i, j - 1))) / dx;
              } else {
                Ey = (this->x(gindex(i, j + 1)) - this->x(gindex(i, j - 1))) /
                     (2 * dx);
              }
              fieldmap(i, j, z) = quadv<double>{0.0, Ex, Ey, Ez};
            }
          }
        }
      });
}

void FieldFDM::save(std::string filename) {
  Logger::GetLogger()->info("FieldFDM: Save matrix in {}", filename);
  std::ofstream file(filename);
  if (file.is_open()) {
    file << x << std::endl;
  }
}

inline int FieldFDM::gindex(int i, int j) const { return i + nx * j; }

} // namespace SpaceCharge
