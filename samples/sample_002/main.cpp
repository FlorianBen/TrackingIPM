#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <random>

#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/field/fields.hpp"
#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/bunch.hpp"

double kernel(double distance, double n) { return exp(-pow(n * distance, 2)); }

int main(int argc, char *argv[]) {
  SpaceCharge::Logger::Init();

  typedef Eigen::Matrix<double, 4, 1> quadv;
  using namespace boost::numeric::odeint;

  if (argc < 2) {
    return 0;
  }

  // Physic constants
  double mproton = 1.6e-27;
  double qelec = 1.6e-19;
  double sol = 3.0e8;

  SpaceCharge::Particle<double> part("proton", 1, SpaceCharge::cst::mproton,
                                     SpaceCharge::cst::lfactor::beta, 0.5);

  SpaceCharge::CSVFileEMField<double> field;
  quadv offset;
  offset << 0.0, 0.0, 0.0, 0.206;
  field.loadEfield(argv[1], offset, 1.0);

  quadv Bfield;
  Bfield << 0.0, 0.0, 0.0, 0.0;

  std::default_random_engine g1(0);
  std::normal_distribution<double> pos_x(0.0, 3e-3);
  std::normal_distribution<double> pos_y(0.0, 3e-3);
  std::uniform_real_distribution<double> pos_z(-0.012, 0.012);
  std::normal_distribution<double> pos_t(-2e-12, 2e-12);

  auto nb_part = 500;
  tbb::concurrent_vector<quadv> pos(nb_part);
  tbb::concurrent_vector<SpaceCharge::state_type2<double>> results_vc;

  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

  for (auto &p : pos) {
    p(0) = 0;
    p(1) = pos_x(g1);
    p(2) = pos_y(g1);
    p(3) = pos_z(g1);
  }

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, nb_part),
      [&](const tbb::blocked_range<size_t> &r) {
        // r particles
        for (auto i = r.begin(); i != r.end(); i++) {
          SpaceCharge::state_type2<double> results;
          // No initial velocity
          quadv vel(.0, .0, .0, .0);
          // Print progress
          if ((i % 1000) == 0)
            std::cout << i << std::endl;
          // Init ODE
          SpaceCharge::state_type2<double> init{pos[i], vel};
          runge_kutta4<SpaceCharge::state_type2<double>> stepper;
          // Define lorentz equation
          auto lorentz = [&](const SpaceCharge::state_type2<double> &x,
                             SpaceCharge::state_type2<double> &dxdt, const double t) {
            quadv Efield;
            Efield << 0.0, 0.0e5, 0.0e5, 0.0;
            // EfieldG.interpolateDelaunay3D_1(x[0], Efield);
            // EfieldG.interpolateNN(x[0], Efield);
            field.interpolateRBF(x[0], Efield, 7, 1, kernel);
            dxdt[0] = x[1];
            dxdt[1] = (qelec / (mproton *
                                sqrt(1 - SpaceCharge::scalar_prod(x[1], x[1]) /
                                             (sol * sol)))) *
                      (Efield + SpaceCharge::vect_prod(x[1], Bfield));
          };
          // Define observer
          auto observer = [&](const SpaceCharge::state_type2<double> &x, double t) {
            if ((abs(x[0][1]) <= 0.05) && (abs(x[0][2]) <= 0.05))
              // Record only if particle is inside IPM
              results.push_back(x[0]);
          };
          // Solve ODE
          integrate_const(stepper, lorentz, init, 0.0, 160e-9, 0.1e-9,
                          observer);
          // Push Track in result vector
          results_vc.push_back(results);
        }
      });
  std::cout << "Tracking done !" << std::endl;
  return 0;
}
