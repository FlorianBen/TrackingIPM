#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <string>

#include "definitions.hpp"

namespace SpaceCharge {

template <class T> class Particle {
private:
  int pcharge;
  T pmass;
  std::string name;
  T pspeed;
  T pec;
  T beta;
  T gamma;

  void updateParticle(cst::lfactor factor, T lfactor_v);
  void updateCGamma();
  void updateCBeta();
  void updateBetaGamma();
  void updateCBetaGamma();

public:
  Particle(std::string name, int p_charge, T pmass, cst::lfactor factor,
           T lfactor_v);

  T getMass() const;
  T getCharge() const;
  T getEc() const;
  T getSpeed() const;
  T getBeta() const;
  T getGamma() const;

  void setEc(const T Ec);
  void setSpeed(const T Speed);
  void setBeta(const T Beta);
  void setGamma(const T Gamma);
};

} // namespace SpaceCharge

#endif