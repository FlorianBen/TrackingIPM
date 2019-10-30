#ifndef BUNCH_HPP
#define BUNCH_HPP

#include <string>

#include "SpaceCharge/definitions.hpp"
#include "SpaceCharge/particle.hpp"

namespace SpaceCharge {

template <class T> class Bunch {
private:
  Particle<T> particle;
  T Qb;
  T ib;
  int pb;
  T rf_freq;
  T dt;

  void check_particle_def(std::string pdg_name);

public:
  Bunch();
  Bunch(std::string pdg_name, T ib, T dt);
  Bunch(Particle<T> particle, T ib, T dt);
  Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib, T dt);

  // virtual T PotentialAt() const = 0;
};

template <class T> class GaussianBunch : public Bunch<T> {};

template <class T> class FEMBunch : public Bunch<T> {};
} // namespace SpaceCharge

#endif