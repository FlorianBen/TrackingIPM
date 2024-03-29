/**
 * \file units.hpp
 * \brief This file contains the definitions of several pysical units.
 * Each function of this project returns SI values, these values can be
 * converted thanks to the definitions contained in this file.
 **/
#ifndef UNITS_HPP
#define UNITS_HPP

#include "definitions.hpp"

/**
 * \namespace SpaceCharge::uni
 * \brief This namespace contains the definitions of pysical units nested in
 * subnamespace.
 **/
namespace SpaceCharge::uni {

const double femto = 1e-15;
const double pico = 1e-12;
const double nano = 1e-9;
const double micro = 1e-6;
const double milli = 1e-3;
const double kilo = 1e3;
const double mega = 1e6;
const double giga = 1e9;
const double tera = 1e12;
const double peta = 1e15;

/**
 * \namespace SpaceCharge::uni::t
 * \brief This namespace contains the definitions of pysical units
 * about time information.
 **/
namespace t {
/**
 * \brief Time base SI unit: second.
 **/
const double s = 1.0;
const double fs = s * femto;
const double ps = s * pico;
const double ns = s * nano;
const double us = s * micro;
const double ms = s * milli;
const double ks = s * kilo;
const double Ms = s * mega;
const double Gs = s * giga;
const double Ts = s * tera;
const double Ps = s * peta;

/**
 * \brief Frequency base unit: Hertz.
 * \f$\frac{1}{T}\f$
 **/
const double Hz = 1.0;
const double fHz = Hz / femto;
const double pHz = Hz / pico;
const double nHz = Hz / nano;
const double uHz = Hz / micro;
const double mHz = Hz / milli;
const double kHz = Hz / kilo;
const double MHz = Hz / mega;
const double GHz = Hz / giga;
const double THz = Hz / tera;
const double PHz = Hz / peta;
} // namespace t

/**
 * \namespace SpaceCharge::uni::d
 * \brief This namespace contains the definitions of pysical units
 * about distance information.
 **/
namespace d {
/**
 * \brief Distance base SI unit: meter.
 **/
const double m = 1.0;
const double fm = m * femto;
const double pm = m * pico;
const double nm = m * nano;
const double um = m * micro;
const double mm = m * milli;
const double km = m * kilo;
const double Mm = m * mega;
const double Gm = m * giga;
const double Tm = m * tera;
const double Pm = m * peta;

} // namespace d

/**
 * \namespace SpaceCharge::uni::m
 * \brief This namespace contains the definitions of pysical units
 * about mass information.
 **/
namespace m {
/**
 * \brief Mass base SI unit: kilogramme.
 **/
const double g = 0.001;
const double fg = g * femto;
const double pg = g * pico;
const double ng = g * nano;
const double ug = g * micro;
const double mg = g * milli;
const double kg = g * kilo;
const double Mg = g * mega;
const double Gg = g * giga;
const double Tg = g * tera;
const double Pg = g * peta;
} // namespace m

/**
 * \namespace SpaceCharge::uni::p
 * \brief This namespace contains the definitions of pysical units
 * about power information.
 **/
namespace p {
/**
 * \brief Power base SI unit: Joule.
 **/
const double J = 1.0;
const double fJ = J * femto;
const double pJ = J * pico;
const double nJ = J * nano;
const double uJ = J * micro;
const double mJ = J * milli;
const double kJ = J * kilo;
const double MJ = J * mega;
const double GJ = J * giga;
const double TJ = J * tera;
const double PJ = J * peta;

/**
 * \brief eV: electronVolt.
 * \f$E_{eV}=E_{Joule} \times 1.602176 \cdot 10^{-19}\f$
 **/
const double eV = 1.0 * cst::e;
const double feV = eV * femto;
const double peV = eV * pico;
const double neV = eV * nano;
const double ueV = eV * micro;
const double meV = eV * milli;
const double keV = eV * kilo;
const double MeV = eV * mega;
const double GeV = eV * giga;
const double TeV = eV * tera;
const double PeV = eV * peta;
} // namespace p

} // namespace SpaceCharge::uni

#endif