#ifndef BUNCH_HPP
#define BUNCH_HPP

#include <string>

#undef Success
#include <Eigen/Dense>

#include <gsl/gsl_integration.h>

#include "SpaceCharge/definitions.hpp"
#include "SpaceCharge/fields.hpp"
#include "SpaceCharge/particle.hpp"
#include "SpaceCharge/point_cloud.hpp"
#include "SpaceCharge/units.hpp"

namespace SpaceCharge {

template <class T> class BunchEMField;

/**
 * \class Bunch bunch.hpp
 * \brief Virutal class that describes a bunch of particles.
 **/
template <class T> class Bunch : public EMField<T> {
  typedef Eigen::Matrix<T, 4, 4> lt_matrix;
  friend BunchEMField<T>;

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
  /**
   * \brief Constructor using an existing particle object.
   * \param[in] particle Particle object.
   * \param[in] ib Particle charge.
   * \param[in] dt Bunch period.
   * \param[in] dir Direction of propagation.
   **/
  Bunch(Particle<T> particle, T ib, T dt, cst::dir dir = cst::dir::z);
  /**
   * \brief Constructor.
   * \param[in] p_charge Particle charge.
   * \param[in] pmass Particle mass.
   * \param[in] p_charge Particle charge.
   * \param[in] factor Lorentz factor type.
   * \param[in] lfactor_v Lorentz value.
   * \param[in] dt Bunch period.
   * \param[in] dir Direction of propagation.
   **/
  Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib, T dt,
        cst::dir dir = cst::dir::z);
  virtual ~Bunch();

  /**
   * \brief Get the particle type used for this bunch.
   **/
  Particle<T> getParticle() const;
  /**
   * \brief Get the bunch current.
   * \return current (\f$A\f$)
   **/
  T getCurrent() const;
  /**
   * \brief Get the number of charge in this bunch.
   * \return charges (\f$C\f$)
   **/
  T getCharges() const;
  /**
   * \brief Get the number of charged particles in this bunch.
   * \return number of particles
   **/
  int getNbParticles() const;
  /**
   * \brief Get the radio frenquency of the bunch.
   * \return number or particle (\f$Hz\f$)
   **/
  T getFreqRF() const;
  /**
   * \brief Get the periodicity of the bunch.
   * \return period (\f$s\f$)
   **/
  T getBunchPeriod() const;

  /**
   * \brief Calculate the Electrical potential at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical potential (\f$V\f$)
   * The derived class must implement this function.
   **/
    virtual T potentialAt(quadv<T> quad) const = 0;
   
  
};

/**
 * \class GaussianBunch bunch.hpp
 * \brief A class that represents a gaussian bunch.
 * The bunch consist of particles gaussianly distributed in 3D.
 * The model description is available
 * <a
 *href="https://www.desy.de/~mpywar/paper/2010/Internal_Report_M_10-01.pdf">here</a>
 * This is an analytic model.
 * The static E field of the bunch is computed in the moving frame, then the EM
 * fields are recovered thank to the Lorentz transform.
 **/
template <class T> class GaussianBunch : public Bunch<T> {
  friend BunchEMField<T>;

private:
  void updateSigma();

  state_type2<T> internalField1() const;
  void internalField2() const;

protected:
  T kurt;
  T sigma0;
  quadv<T> sigma;
  quadv<T> sigma_;

public:
  /**
   * \brief Constructor using an existing particle object.
   * \param[in] particle Particle object.
   * \param[in] ib Particle charge.
   * \param[in] dt Bunch period.
   * \param[in] dir Direction of propagation.
   **/
  GaussianBunch(Particle<T> particle, T ib, T dt, cst::dir dir = cst::dir::z);
  /**
   * \brief Constructor.
   * \param[in] p_charge Particle charge.
   * \param[in] pmass Particle mass.
   * \param[in] p_charge Particle charge.
   * \param[in] factor Lorentz factor type.
   * \param[in] lfactor_v Lorentz value.
   * \param[in] dt Bunch period.
   * \param[in] dir Direction of propagation.
   **/
  GaussianBunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt, cst::dir dir = cst::dir::z);

  /**
   * \brief Get the sigma of the bunch in the rest frame.
   * \return sigma vector (\f$m\f$)
   **/
  quadv<T> getSigma() const;
  /**
   * \brief Get the sigma of the bunch in the moving frame.
   * \return sigma vector (\f$m\f$)
   **/
  quadv<T> getSigmaB() const;
  /**
   * \brief Set the sigma of the bunch in the rest frame.
   * \param[in] sigma sigma quadvector (\f$m\f$)
   **/
  void setSigma(quadv<T> sigma);
  /**
   * \brief Set the sigma of the bunch in the rest frame.
   * \param[in] sigma sigma vector (\f$m\f$)
   **/
  void setSigma(Eigen::Matrix<T, 3, 1> sigma);
  /**
   * \brief Set the sigma of the bunch in the rest frame.
   * \param[in] index direction
   * \param[in] sigma sigma value (\f$m\f$)
   **/
  void setSigma(int index, T sigma);
  /**
   * \brief Set the position of the bunch in the rest frame.
   * \param[in] sigma position quadvector (\f$m\f$)
   **/
  //  void setPosition(Eigen::Matrix<T, 4, 1> pos);
  /**
   * \brief Set the position of the bunch in the rest frame.
   * \param[in] sigma position vector (\f$m\f$)
   **/
  //  void setPosition(Eigen::Matrix<T, 3, 1> pos);
  /**
   * \brief Set the position of the bunch in the rest frame.
   * \param[in] index direction
   * \param[in] sigma position value (\f$m\f$)
   **/
  //  void setPosition(int index, T pos);
  /**
   * \brief Calculate the Electrical potential at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical potential (\f$V\f$)
   **/
  T potentialAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector (\f$V/m\f$)
   **/
  quadv<T> EfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Magnetic field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector (\f$T\f$)
   **/
  quadv<T> MagfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   * The derived class must implement this function.
   **/
  virtual state_type2<T> EMfielddAt(quadv<T> quad) const override;

  state_type2<T> transformFields(state_type2<T> fields) const;
};

template <class T> class FEMBunch : public Bunch<T> {};

/**
 * \class gsl_function_pp bunch.hpp
 * \brief Helper class that wraps a lambda/functor into a GSL function.
 * This class helps to integrate with GSL directly with a functor or a lambda
 * expression.
 **/
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

/**
 * \class FieldBunch fields.hpp
 * \brief A class that represents an EM field created by particle bunches.
 * Several bunch can be added. The resulting EM field is the sum of the
 * contribution of each bunch.
 **/
template <class T> class BunchEMField : public EMField<T> {
private:
  std::vector<std::unique_ptr<Bunch<T>>> bunches;
  bool use_periodicity;
  T local_time;

public:
  BunchEMField();

  /**
   * \brief Add a bunch as source of EM field.
   * \param[in] bunch Pointer to a bunch.
   **/
  void addBunch(std::unique_ptr<Bunch<T>> bunch);
  /**
   * \brief Use the period of the bunch.
   * \param[in] use Use periodicity if True.
   **/
  void usePeriodicity(bool use = true);
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector.
   **/
  virtual quadv<T> EfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector.
   **/
  virtual quadv<T> MagfieldAt(quadv<T> quad) const override;

  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   * The derived class must implement this function.
   **/
  virtual state_type2<T> EMfielddAt(quadv<T> quad) const override;
};

} // namespace SpaceCharge

#endif