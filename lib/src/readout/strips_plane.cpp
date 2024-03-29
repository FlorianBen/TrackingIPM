#include "SpaceCharge/readout/strips_plane.hpp"
#include "SpaceCharge/core/alogger.hpp"

#include <fstream>
#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>

// TODO: Implement dim gap independant version.
// TODO: Add a way to select Eigen solver.
// TODO: Add a templated version.
// TODO: Better documentation.

namespace SpaceCharge {

// PUBLIC METHOD

StripsPlane::StripsPlane(const std::string filename) : gap(0.01) {
  readFile(filename);
  computeBestSize();
}

StripsPlane::StripsPlane(const std::string filename,
                         const SpaceCharge::quadv<size_t> sizes)
    : gap(0.01) {
  readFile(filename);
  dx = pcb_size / sizes(1);
  dy = pcb_size / sizes(1);

  nx = sizes(1);
  ny = sizes(2);

  b.resize(nx * ny);
  mat.resize(nx * ny, nx * ny);
}

void StripsPlane::solvePotential(const int strip_on) {
  initMatrix(strip_on);
  solve();
}

void StripsPlane::getPotential(FieldMap<double> &fieldmap, int zindex) const {
  tbb::parallel_for(
      tbb::blocked_range3d<int>(0, 1, 0, ny, 0, nx),
      [&](const tbb::blocked_range3d<int> &r) {
        for (int y = r.rows().begin(), y_end = r.rows().end(); y < y_end; y++) {
          for (int x = r.cols().begin(), x_end = r.cols().end(); x < x_end;
               x++) {
            for (int z = r.pages().begin(), z_end = r.pages().end(); z < z_end;
                 z++) {
              fieldmap(x, y, zindex) =
                  quadv<double>{0.0, this->x(x + nx * y), 0.0, 0.0};
            }
          }
        }
      });
}

void StripsPlane::getField(FieldMap<double> &fieldmap, int zindex) const {
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
              fieldmap(i, j, zindex) = quadv<double>{0.0, Ex, Ey, Ez};
            }
          }
        }
      });
}

size_t StripsPlane::getSizeX() const { return nx; }

size_t StripsPlane::getSizeY() const { return ny; }

double StripsPlane::getGapX() const { return dx; }

double StripsPlane::getGapY() const { return dy; }

size_t StripsPlane::getNbStrips() const { return total_strips; }

// PRIVATE METHODS

void StripsPlane::readFile(const std::string filename) {
  std::ifstream infile(filename);
  double a, b;

  while (infile >> a >> b) {
    a *= 1e-3;
    b *= 1e-3;
    std::pair<double, double> pair_corr{a, b};
    strips_pairs.push_back(pair_corr);
  }

  auto x_min = std::begin(strips_pairs)->first;
  auto x_max = std::prev(std::end(strips_pairs))->second;

  pcb_size = x_max - x_min;

  auto min_inter = pcb_size;
  auto min_strips = pcb_size;
  auto min_dist = pcb_size;
  auto old_first = 0.0;

  for (auto &pair : strips_pairs) {
    pair.first = pair.first - x_min;
    pair.second = pair.second - x_min;

    auto min_intert = pair.second - old_first;
    if (min_intert < min_inter) {
      min_inter = min_intert;
    }

    auto min_stripst = pair.second - pair.first;
    if (min_stripst < min_strips) {
      min_strips = min_stripst;
    }

    if (min_strips < min_inter) {
      min_dx = min_strips;
    } else {
      min_dx = min_inter;
    }

    old_first = pair.first;
  }
  total_strips = strips_pairs.size() - 2;
}

void StripsPlane::initMatrix(const int strip_on) {
  std::vector<Tri> coefficients_t;
  auto strip_nb = 0;
  b.setZero();
  mat.setZero();
  mat.data().squeeze();
  for (auto j = 0; j < ny; j++) {
    for (auto i = 0; i < nx; i++) {
      auto indl = gindex(i, j);
      if (j == 0) {
        b(indl) = -0;
        coefficients_t.push_back(Tri(indl, indl, 1));
        continue;
      }
      if (j == (ny - 1)) {
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
          coefficients_t.push_back(Tri(indl, gindex(i, j - 1), 1));
          coefficients_t.push_back(Tri(indl, indl, -2));
        }
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
  mat.setFromTriplets(coefficients_t.begin(), coefficients_t.end());
}

void StripsPlane::solve() {
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(1e-7);
  solver.compute(mat);
  x = solver.solve(b);
  Logger::GetLogger()->info("FieldFDM: Solving done in {} iterations, error {}",
                            solver.iterations(), solver.error());
}

bool StripsPlane::isStrips(const int i, int &stripnb) const {
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

void StripsPlane::computeBestSize() {
  auto sizet = 200;

  while ((min_dx / 20) < (pcb_size / sizet)) {
    sizet += 500;
  }

  dx = pcb_size / sizet;
  dy = pcb_size / sizet;

  nx = sizet;
  ny = gap / dy;

  b.resize(nx * ny);
  mat.resize(nx * ny, nx * ny);
}

inline int StripsPlane::gindex(int i, int j) const { return i + nx * j; }

} // namespace SpaceCharge
