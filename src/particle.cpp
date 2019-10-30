#include "SpaceCharge/particle.hpp"

namespace SpaceCharge {

// template <class T> Particle<T>::Particle(){};

template <class T>
Particle<T>::Particle(std::string name, int p_charge, T pmass,
                      cst::lfactor factor, T lfactor_v)
    : name(name), pcharge(p_charge), pmass(pmass) {
  updateParticle(factor, lfactor_v);
};

template <class T>
void Particle<T>::updateParticle(cst::lfactor factor, T lfactor_v) {
  switch (factor) {
  case cst::lfactor::beta:
    beta = lfactor_v;
    updateCGamma();
    break;
  case cst::lfactor::gamma:
    gamma = lfactor_v;
    updateCBeta();
    break;
  case cst::lfactor::speed:
    pspeed = lfactor_v;
    updateBetaGamma();
    break;
  case cst::lfactor::ec:
    pec = lfactor_v;
    updateCBetaGamma();
    break;
  default:
    break;
  }
};

template <class T> void Particle<T>::updateCGamma() {
  pspeed = beta * cst::sol;
  gamma = 1.0 / std::sqrt(1.0 - (std::pow(beta, 2)));
  pec = (1 - gamma) * pmass * std::pow(cst::sol, 2);
};

template <class T> void Particle<T>::updateCBeta() {
  pec = (1 - gamma) * pmass * std::pow(cst::sol, 2);
  beta = std::sqrt(1.0 - (1.0 / std::pow(gamma, 2)));
  pspeed = cst::sol * beta;
};

template <class T> void Particle<T>::updateBetaGamma() {
  beta = pspeed / cst::sol;
  gamma = 1.0 / std::sqrt(1.0 - (std::pow(beta, 2)));
  pec = (1 - gamma) * pmass * std::pow(cst::sol, 2);
};

template <class T> void Particle<T>::updateCBetaGamma() {
  gamma = pec / (pmass * std::pow(cst::sol, 2)) + 1;
  beta = std::sqrt(1.0 - (1.0 / std::pow(gamma, 2)));
  pspeed = beta * cst::sol;
};

template <class T> T Particle<T>::getMass() const { return pmass; };

template <class T> T Particle<T>::getCharge() const { return pcharge; };

template <class T> T Particle<T>::getEc() const { return pec; };

template <class T> T Particle<T>::getSpeed() const { return pspeed; };

template <class T> T Particle<T>::getBeta() const { return beta; };

template <class T> T Particle<T>::getGamma() const { return gamma; };

template <class T> void Particle<T>::setEc(const T Ec) {
  pec = Ec;
  updateCBetaGamma();
};

template <class T> void Particle<T>::setSpeed(const T Speed) {
  pspeed = Speed;
  updateBetaGamma();
};

template <class T> void Particle<T>::setBeta(const T Beta) {
  beta = Beta;
  updateCGamma();
};

template <class T> void Particle<T>::setGamma(const T Gamma) {
  gamma = Gamma;
  updateCBeta();
};

template class Particle<double>;
template class Particle<float>;

} // namespace SpaceCharge