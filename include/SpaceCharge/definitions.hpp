#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <boost/math/constants/constants.hpp>

namespace SpaceCharge {
/**
 * \namespace SpaceCharge::cst
 * \brief This namespace contains the definitions of some pysical constants.
 **/
namespace cst {
  // TODO: Put unit for Doxygen
//! Definition of \f$\pi\f$ as constexpr
constexpr double pi = boost::math::constants::pi<double>();
//! Gravitational acceleration
const double g = 9.8;
//! Boltzmann constant
const double k = 1.380649e-23;
//! Avogadro constant
const double Na = 6.02214076e23;
//! Perfect law constant
const double R = Na * k;
//! Usual space directions
enum dir { x, y, z };
//! Usual ways fro defining the Lorentz factors
enum lfactor { speed, beta, gamma, ec };
//! Plank constant
const double hp = 6.62607015e-34;
//! Reduced Plank constant
const double hpb = hp / (2 * pi);
//! Speed of light (vacuum) in SI
const double sol = 2.99792458e8;
//! Elementary charge in SI
const double e = 1.6021766208e-19;
//! Vacuum permeability
const double mu0 = 4 * pi * 1e-7;
//! Vacuum permittivity
const double eps0 = 1 / (mu0 * sol * sol);

//! Mass of one proton in SI
const double mproton = 1.672621898e-27;
//! Mass of one electron in SI
const double melectron = 9.10938356e-31;
//! Mass of one neutron in SI
const double mneutron = 1.674927471e-27;

} // namespace cst
} // namespace SpaceCharge
#endif