#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <string>

#include "definitions.hpp"

namespace SpaceCharge {

/**
 * \class Particle particle.hpp
 * \brief This class describes the particles that compose the bunch.
 *
 **/
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
  /**
   * \brief Constructor.
   * \param[in] name Particle name.
   * \param[in] pcharge Particle charge.
   * \param[in] pmass Particle mass.
   * \param[in] factor Select speed, \f$\beta\f$, \f$\gamma\f$ or \f$E_{c}\f$
   * factor.
   * \param[in] lfactor_v Value of the factor.
   **/
  Particle(std::string name, int pcharge, T pmass, cst::lfactor factor,
           T lfactor_v);

  /**
   * \brief Get Particle mass.
   * \return The value of the particle mass in kg.
   **/
  T getMass() const;
  /**
   * \brief Get Particle charge.
   * \return The value of the particle charge in C.
   **/
  T getCharge() const;
  /**
   * \brief Get Particle kinetic energy.
   * \return The value of the particle mass in J.
   **/
  T getEc() const;
  /**
   * \brief Get Particle speed.
   * \return The value of the particle mass in m/s.
   **/
  T getSpeed() const;
  /**
   * \brief Get Particle \f$\beta\f$.
   * \return The value of the particle \f$\beta\f$.
   **/
  T getBeta() const;
  /**
   * \brief Get Particle \f$\gamma\f$.
   * \return The value of the particle \f$\gamma\f$.
   **/
  T getGamma() const;

  /**
   * \brief
   * \param[in]
   **/
  void setEc(const T Ec);
  /**
   * \brief
   * \param[in]
   **/
  void setSpeed(const T Speed);
  /**
   * \brief
   * \param[in]
   **/
  void setBeta(const T Beta);
  /**
   * \brief
   * \param[in]
   **/
  void setGamma(const T Gamma);
};

} // namespace SpaceCharge

#endif