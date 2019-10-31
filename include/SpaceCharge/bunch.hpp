#ifndef BUNCH_HPP
#define BUNCH_HPP

#include <string>

#include <Eigen/Dense>

#include "SpaceCharge/definitions.hpp"
#include "SpaceCharge/particle.hpp"

namespace SpaceCharge {

template <class T> class Bunch {
protected:
  Particle<T> particle;
  T Qb;
  T ib;
  int pb;
  T rf_freq;
  T dt;

  void check_particle_def(std::string pdg_name);
  void updateBunch();

public:
  // Bunch(std::string pdg_name, T ib, T dt);
  Bunch(Particle<T> particle, T ib, T dt);
  Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib, T dt);

  Particle<T> getParticle() const;
  T getCurrent() const;
  T getCharges() const;
  int getNbParticles() const;
  T getFreqRF() const;
  T getBunchPeriod() const;

  // virtual T PotentialAt() const = 0;
};

template <class T> class GaussianBunch : public Bunch<T> {
  typedef Eigen::Matrix<T, 4, 1> quadv;

protected:
  quadv sigma;
  quadv sigma_;

public:
  GaussianBunch(Particle<T> particle, T ib, T dt, quadv sigma);
  GaussianBunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt, quadv sigma);
};

template <class T> class FEMBunch : public Bunch<T> {};
} // namespace SpaceCharge

#endif