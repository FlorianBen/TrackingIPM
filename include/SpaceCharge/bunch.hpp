#ifndef BUNCH_HPP
#define BUNCH_HPP

#include <string>

#undef Success
#include <Eigen/Dense>

#include <gsl/gsl_integration.h>

#include "SpaceCharge/definitions.hpp"
#include "SpaceCharge/particle.hpp"
#include "SpaceCharge/units.hpp"

namespace SpaceCharge {

template <class T> class FieldBunch;

template <class T> class Bunch {
  typedef Eigen::Matrix<T, 4, 1> quadv;
  typedef Eigen::Matrix<T, 4, 4> lt_matrix;
  friend FieldBunch<T>;

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
  virtual ~Bunch();

  Particle<T> getParticle() const;
  T getCurrent() const;
  T getCharges() const;
  int getNbParticles() const;
  T getFreqRF() const;
  T getBunchPeriod() const;

  virtual T potentialAt(quadv quad) = 0;
  virtual quadv EfieldAt(quadv quad) = 0;
  virtual quadv MagfieldAt(quadv quad) = 0;
};

template <class T> class GaussianBunch : public Bunch<T> {
  typedef Eigen::Matrix<T, 4, 1> quadv;
  friend FieldBunch<T>;

private:
  void updateSigma();
  void updatePosition();
  void updateFields();
  void internalPotential();
  void internalField1();
  void internalField2();

protected:
  T kurt;
  T sigma0;
  quadv sigma;
  quadv sigma_;
  quadv pos;
  quadv pos_;
  T V;
  T error_V;
  quadv E;
  quadv E_;
  quadv B;
  quadv B_;

public:
  GaussianBunch(Particle<T> particle, T ib, T dt, cst::dir dir = cst::dir::z);
  GaussianBunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt, cst::dir dir = cst::dir::z);

  Eigen::Matrix<T, 4, 1> getSigma() const;
  Eigen::Matrix<T, 4, 1> getSigmaB() const;

  void setSigma(Eigen::Matrix<T, 4, 1> sigma);
  void setSigma(Eigen::Matrix<T, 3, 1> sigma);
  void setSigma(int index, T sigma);

  void setPosition(Eigen::Matrix<T, 4, 1> pos);
  void setPosition(Eigen::Matrix<T, 3, 1> pos);
  void setPosition(int index, T pos);

  T potentialAt(quadv quad) override;
  quadv EfieldAt(quadv quad) override;
  quadv MagfieldAt(quadv quad) override;
};

template <class T> class FEMBunch : public Bunch<T> {};

template <typename F> class gsl_function_pp : public gsl_function {
public:
  gsl_function_pp(const F &func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params = this;
  }

private:
  const F &_func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp *>(params)->_func(x);
  }
};

} // namespace SpaceCharge

#endif