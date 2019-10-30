#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <string>

#include "SpaceCharge/definitions.hpp"

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

public:
  Particle(std::string name, int p_charge, T pmass, cst::lfactor factor,
           T lfactor_v);
};

} // namespace SpaceCharge

#endif