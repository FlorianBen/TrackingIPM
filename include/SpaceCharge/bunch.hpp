#ifndef BUNCH_HPP
#define BUNCH_HPP

#include <string>

#include <Eigen/Dense>

#include "SpaceCharge/definitions.hpp"
#include "SpaceCharge/particle.hpp"
#include "SpaceCharge/units.hpp"

namespace SpaceCharge {

template <class T> class Bunch {
  typedef Eigen::Matrix<T, 4, 4> lt_matrix;

protected:
  lt_matrix lboost;
  lt_matrix ltransform;
  cst::dir dir;
  Particle<T> particle;
  T Qb;
  T ib;
  int pb;
  T rf_freq;
  T dt;

  void check_particle_def(std::string pdg_name);
  void updateBunch();
  void updateBoost();

public:
  // Bunch(std::string pdg_name, T ib, T dt);
  Bunch(Particle<T> particle, T ib, T dt, cst::dir dir = cst::dir::z);
  Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib, T dt,
        cst::dir dir = cst::dir::z);

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

private:
  void updateSigma();

protected:
  quadv sigma;
  quadv sigma_;
  quadv pos;
  quadv pos_;

public:
  GaussianBunch(Particle<T> particle, T ib, T dt, cst::dir dir = cst::dir::z);
  GaussianBunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt, cst::dir dir = cst::dir::z);

  Eigen::Matrix<T, 4, 1> getSigma() const;
  Eigen::Matrix<T, 4, 1> getSigmaB() const;

  void setSigma(Eigen::Matrix<T, 4, 1> sigma);
  void setSigma(Eigen::Matrix<T, 3, 1> sigma);
  void setSigma(int index, T sigma);
};

template <class T> class FEMBunch : public Bunch<T> {};
} // namespace SpaceCharge

#endif