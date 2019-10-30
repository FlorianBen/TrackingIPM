#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <boost/math/constants/constants.hpp>

namespace SpaceCharge::cst {
constexpr double pi = boost::math::constants::pi<double>();

const double g = 9.8;

const double k = 1.380649e-23;
const double Na = 6.02214076e23;
const double R = Na * k;

enum lfactor { speed, beta, gamma };
const double hp = 6.62607015e-34;
const double hpb = hp / (2 * pi);
const double sol = 2.99792458e8;
const double e = 1.6021766208e-19;
const double mu0 = 4 * pi * 1e-7;
const double eps0 = 1 / (mu0 * sol * sol);

//TODO: Add some constant mass value 
const double mproton =0;

} // namespace SpaceCharge::cst

#endif