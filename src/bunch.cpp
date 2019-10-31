#include "SpaceCharge/bunch.hpp"

namespace SpaceCharge {

// template <class T>
// Bunch<T>::Bunch(std::string pdg_name, T ib, T dt) : ib(ib), dt(dt) {
//   check_particle_def(pdg_name);
// };

template <class T>
Bunch<T>::Bunch(Particle<T> particle, T ib, T dt)
    : particle(particle), ib(ib), dt(dt) {
  updateBunch();
};

template <class T>
Bunch<T>::Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt)
    : particle("particle", p_charge, pmass, factor, lfactor_v), ib(ib), dt(dt) {
  updateBunch();
};

template <class T>
void Bunch<T>::check_particle_def(std::string pdg_name){

};

template <class T> void Bunch<T>::updateBunch() {
  Qb = ib * dt;
  pb = Qb / cst::e;
  rf_freq = 1.0 / dt;
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
GaussianBunch<T>::GaussianBunch(Particle<T> particle, T ib, T dt, quadv sigma)
    : Bunch<T>(particle, ib, dt), sigma(sigma){};

template <class T>
GaussianBunch<T>::GaussianBunch(int p_charge, T pmass, cst::lfactor factor,
                                T lfactor_v, T ib, T dt, quadv sigma)
    : Bunch<T>(p_charge, pmass, factor, lfactor_v), sigma(sigma){};

template class Bunch<double>;
template class Bunch<float>;

} // namespace SpaceCharge