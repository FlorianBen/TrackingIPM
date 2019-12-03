#include "SpaceCharge/bunch.hpp"

#include <iostream>

namespace SpaceCharge {

// template <class T>
// Bunch<T>::Bunch(std::string pdg_name, T ib, T dt) : ib(ib), dt(dt) {
//   check_particle_def(pdg_name);
// };

template <class T>
Bunch<T>::Bunch(Particle<T> particle, T ib, T dt, cst::dir dir)
    : particle(particle), ib(ib), dt(dt), dir(dir) {
  updateBunch();
  updateBoost();
};

template <class T>
Bunch<T>::Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt, cst::dir dir)
    : particle("particle", p_charge, pmass, factor, lfactor_v), ib(ib), dt(dt),
      dir(dir) {
  updateBunch();
  updateBoost();
};

template <class T>
void Bunch<T>::check_particle_def(std::string pdg_name){

};

template <class T> void Bunch<T>::updateBunch() {
  Qb = ib * dt;
  pb = Qb / cst::e;
  rf_freq = 1.0 / dt;
};

template <class T> void Bunch<T>::updateBoost() {
  lboost = Eigen::Matrix<T, 4, 4>::Identity();
  ltransform = Eigen::Matrix<T, 4, 4>::Identity();
  switch (dir) {
  case cst::dir::x:
    lboost(0, 0) = particle.getGamma();
    lboost(1, 1) = particle.getGamma();
    lboost(0, 1) = -particle.getGamma() * particle.getBeta();
    lboost(1, 0) = -particle.getGamma() * particle.getBeta();
    ltransform(1, 1) = particle.getGamma();
    break;
  case cst::dir::y:
    lboost(0, 0) = particle.getGamma();
    lboost(2, 2) = particle.getGamma();
    lboost(0, 2) = -particle.getGamma() * particle.getBeta();
    lboost(2, 0) = -particle.getGamma() * particle.getBeta();
    ltransform(2, 2) = particle.getGamma();
    break;
  case cst::dir::z:
    lboost(0, 0) = particle.getGamma();
    lboost(3, 3) = particle.getGamma();
    lboost(0, 3) = -particle.getGamma() * particle.getBeta();
    lboost(3, 0) = -particle.getGamma() * particle.getBeta();
    ltransform(3, 3) = particle.getGamma();
    break;
  default:
    break;
  }
};

template <class T> Particle<T> Bunch<T>::getParticle() const {
  return Particle<T>(particle);
};

template <class T> T Bunch<T>::getCurrent() const { return ib; };

template <class T> T Bunch<T>::getCharges() const { return Qb; };

template <class T> int Bunch<T>::getNbParticles() const { return pb; };

template <class T> T Bunch<T>::getFreqRF() const { return rf_freq; };

template <class T> T Bunch<T>::getBunchPeriod() const { return dt; };

template <class T>
GaussianBunch<T>::GaussianBunch(Particle<T> particle, T ib, T dt, cst::dir dir)
    : Bunch<T>(particle, ib, dt, dir) {
  sigma(0) = 0.;
  sigma(1) = 3.;
  sigma(2) = 3.;
  sigma(3) = 1.5;
  updateSigma();
};

template <class T>
GaussianBunch<T>::GaussianBunch(int p_charge, T pmass, cst::lfactor factor,
                                T lfactor_v, T ib, T dt, cst::dir dir)
    : Bunch<T>(p_charge, pmass, factor, lfactor_v, ib, dt, dir) {
  updateSigma();
};

template <class T>
Eigen::Matrix<T, 4, 1>  GaussianBunch<T>::getSigma() const{
  return sigma;
};

template <class T>
Eigen::Matrix<T, 4, 1>  GaussianBunch<T>::getSigmaB() const{
  return sigma_;
};

template <class T>
void GaussianBunch<T>::setSigma(Eigen::Matrix<T, 4, 1> sigma) {
  this->sigma = sigma;
  updateSigma();
}

template <class T>
void GaussianBunch<T>::setSigma(Eigen::Matrix<T, 3, 1> sigma) {
  this->sigma.tail(3) = sigma;
  updateSigma();
}

template <class T> void GaussianBunch<T>::setSigma(int index, T sigma) {
  if ((index < 1) || (index > 3)) {
    return;
  }
  this->sigma(index) = sigma;
  updateSigma();
}

template <class T> void GaussianBunch<T>::updateSigma() {
  this->sigma_ = (this->sigma).transpose() * (this->ltransform);
};

template class Bunch<double>;
template class Bunch<float>;

template class GaussianBunch<double>;
template class GaussianBunch<float>;

} // namespace SpaceCharge