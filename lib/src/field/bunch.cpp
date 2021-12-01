#include <spdlog/fmt/ostr.h>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/field/bunch.hpp"


// TODO: Add sigma parameter to the GaussianBunch constructor.

namespace SpaceCharge {

// template <class T>
// Bunch<T>::Bunch(std::string pdg_name, T ib, T dt) : ib(ib), dt(dt) {
//   check_particle_def(pdg_name);
// };

template <class T>
Bunch<T>::Bunch(Particle<T> particle, T ib, T dt, cst::dir dir)
    : particle(particle), ib(ib), dt(dt), dir(dir) {
  updateBunch();
  updateBoost();
};

template <class T>
Bunch<T>::Bunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib,
                T dt, cst::dir dir)
    : particle("particle", p_charge, pmass, factor, lfactor_v), ib(ib), dt(dt),
      dir(dir) {
  updateBunch();
  updateBoost();
};

template <class T> Bunch<T>::~Bunch() {}

template <class T>
void Bunch<T>::check_particle_def(std::string pdg_name){

};

template <class T> void Bunch<T>::updateBunch() {
  Qb = ib * dt;
  pb = Qb / cst::e;
  rf_freq = 1.0 / dt;
};

template <class T> void Bunch<T>::updateBoost() {
  lboost = Eigen::Matrix<T, 4, 4>::Identity();
  ltransform = Eigen::Matrix<T, 4, 4>::Identity();
  switch (dir) {
  case cst::dir::x:
    lboost(0, 0) = particle.getGamma();
    lboost(1, 1) = particle.getGamma();
    lboost(0, 1) = -particle.getGamma() * particle.getBeta();
    lboost(1, 0) = -particle.getGamma() * particle.getBeta();
    ltransform(1, 1) = particle.getGamma();
    break;
  case cst::dir::y:
    lboost(0, 0) = particle.getGamma();
    lboost(2, 2) = particle.getGamma();
    lboost(0, 2) = -particle.getGamma() * particle.getBeta();
    lboost(2, 0) = -particle.getGamma() * particle.getBeta();
    ltransform(2, 2) = particle.getGamma();
    break;
  case cst::dir::z:
    lboost(0, 0) = particle.getGamma();
    lboost(3, 3) = particle.getGamma();
    lboost(0, 3) = -particle.getGamma() * particle.getBeta();
    lboost(3, 0) = -particle.getGamma() * particle.getBeta();
    ltransform(3, 3) = particle.getGamma();
    break;
  default:
    break;
  }
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
std::ostream &operator<<(std::ostream &os, const Bunch<T> &c) {
  return os << "Bunch";
}

template <class T>
GaussianBunch<T>::GaussianBunch(Particle<T> particle, T ib, T dt, cst::dir dir)
    : Bunch<T>(particle, ib, dt, dir) {
  sigma(0) = 0.e-3;
  sigma(1) = 3.e-3;
  sigma(2) = 2.e-3;
  sigma(3) = 1.5e-3;
  updateSigma();
};

template <class T>
GaussianBunch<T>::GaussianBunch(int p_charge, T pmass, cst::lfactor factor,
                                T lfactor_v, T ib, T dt, cst::dir dir)
    : Bunch<T>(p_charge, pmass, factor, lfactor_v, ib, dt, dir) {
  updateSigma();
};

template <class T> Eigen::Matrix<T, 4, 1> GaussianBunch<T>::getSigma() const {
  return sigma;
};

template <class T> Eigen::Matrix<T, 4, 1> GaussianBunch<T>::getSigmaB() const {
  return sigma_;
};

template <class T>
void GaussianBunch<T>::setSigma(Eigen::Matrix<T, 4, 1> sigma) {
  this->sigma = sigma;
  updateSigma();
}

template <class T>
void GaussianBunch<T>::setSigma(Eigen::Matrix<T, 3, 1> sigma) {
  this->sigma.tail(3) = sigma;
  updateSigma();
}

template <class T> void GaussianBunch<T>::setSigma(int index, T sigma) {
  if ((index < 1) || (index > 3)) {
    return;
  }
  this->sigma(index) = sigma;
  updateSigma();
}

template <class T> void GaussianBunch<T>::updateSigma() {
  this->sigma_ = (this->sigma).transpose() * (this->ltransform);
  switch (this->dir) {
  case cst::dir::x:
    sigma0 = std::sqrt(2 * (std::pow(sigma_(3), 2) - std::pow(sigma_(2), 2)));
    kurt = sigma_(2) / sigma_(3);
    break;
  case cst::dir::y:
    sigma0 = std::sqrt(2 * (std::pow(sigma_(1), 2) - std::pow(sigma_(3), 2)));
    kurt = sigma_(3) / sigma_(1);
    break;
  case cst::dir::z:
    sigma0 = std::sqrt(2 * (std::pow(sigma_(1), 2) - std::pow(sigma_(2), 2)));
    kurt = sigma_(2) / sigma_(1);
    break;
  default:
    break;
  }
};

template <class T> T GaussianBunch<T>::potentialAt(quadv<T> quad) const {
  T V;
  quadv<T> pos = quad;
  quadv<T> pos_ = pos.transpose() * this->lboost;

  auto factor = this->Qb / (4.0 * cst::pi * cst::eps0 * std::sqrt(cst::pi));
  auto fun = [&](T q) -> T {
    Eigen::Matrix<T, 4, 1> qv;
    qv << 0.0, q + 2 * std::pow(sigma_(1), 2), q + 2 * std::pow(sigma_(2), 2),
        q + 2 * std::pow(sigma_(3), 2);
    auto term_root = std::sqrt(qv(1) * qv(2) * qv(3));
    auto term_rat = -((std::pow(pos_(1), 2)) / qv(1)) -
                    ((std::pow(pos_(2), 2)) / qv(2)) -
                    ((std::pow(pos_(3), 2)) / qv(3));
    return std::exp(term_rat) / term_root;
  };
  gsl_function_pp<decltype(fun)> Fp(fun);
  gsl_function *F = static_cast<gsl_function *>(&Fp);
  gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(1000);

  T abs_error = 1.0e-8;
  T rel_error = 1.0e-8;
  T result;
  T error;

  gsl_integration_qagiu(F, 0.0, abs_error, rel_error, 1000, work_ptr, &result,
                        &error);
  V = result;
  auto error_V = error;
  gsl_integration_workspace_free(work_ptr);

  return V;
};

template <class T> quadv<T> GaussianBunch<T>::EfieldAt(quadv<T> quad) const {
  return EMfieldAt(quad)[0];
};

template <class T> quadv<T> GaussianBunch<T>::MagfieldAt(quadv<T> quad) const {
  return EMfieldAt(quad)[1];
};

template <class T>
state_type2<T> GaussianBunch<T>::EMfieldAt(quadv<T> quad) const {
  quadv<T> E_{0.0, 0.0, 0.0, 0.0};
  quadv<T> B_{0.0, 0.0, 0.0, 0.0};
  quadv<T> pos = quad;
  quadv<T> pos_ = pos.transpose() * this->lboost;

  auto factor = this->Qb / (2.0 * cst::eps0 * std::sqrt(cst::pi) * sigma0) * this->particle.getCharge();
  auto axis = 0;

  auto fun = [&](T eps) -> T {
    Eigen::Matrix<T, 3, 3> factor_int;
    Eigen::Matrix<T, 3, 3> qe;
    qe << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        std::pow(sigma0, 2) * (1.0 / (1.0 - std::pow(eps, 2))),
        std::pow(sigma0, 2) * (std::pow(eps, 2) / (1.0 - std::pow(eps, 2))),
        2 * (std::pow(sigma_(this->dir + 1), 2) - std::pow(sigma_(1), 2)) +
            std::pow(sigma0, 2) * (1.0 / (1.0 - std::pow(eps, 2)));
    factor_int << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 / std::pow(eps, 2),
        qe(this->dir, 0) / qe(this->dir, 2);
    auto term_root = std::sqrt(qe(2, this->dir));
    auto term_rat = -((std::pow(pos_(1), 2)) / qe(this->dir, 0)) -
                    ((std::pow(pos_(2), 2)) / qe(this->dir, 1)) -
                    ((std::pow(pos_(3), 2)) / qe(this->dir, 2));
    return factor_int(this->dir, axis) * std::exp(term_rat) / term_root;
  };

  gsl_function_pp<decltype(fun)> Fp(fun);
  gsl_function *F = static_cast<gsl_function *>(&Fp);
  gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(1000);
  T abs_error = 1.0e-8;
  T rel_error = 1.0e-8;
  T result;
  T error;

  for (auto i = 0; i < 3; i++) {
    axis = i;
    gsl_integration_qag(F, kurt, 1.0, abs_error, rel_error, 1000,
                        GSL_INTEG_GAUSS31, work_ptr, &result, &error);
    E_(i + 1) = 2.0 * factor * (pos_(i + 1) / (sigma0 * cst::pi)) * result;
  }
  gsl_integration_workspace_free(work_ptr);
  return transformFields(state_type2<T>{E_, B_});
}

template <class T> state_type2<T> GaussianBunch<T>::internalField1() const {
  // auto factor = this->Qb / (2.0 * cst::eps0 * std::sqrt(cst::pi) * sigma0);
  // auto axis = 0;
  // auto fun = [&](T eps) -> T {
  //   Eigen::Matrix<T, 3, 3> factor_int;
  //   Eigen::Matrix<T, 3, 3> qe;

  //   qe << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  //       std::pow(sigma0, 2) * (1.0 / (1.0 - std::pow(eps, 2))),
  //       std::pow(sigma0, 2) * (std::pow(eps, 2) / (1.0 - std::pow(eps, 2))),
  //       2 * (std::pow(sigma_(this->dir + 1), 2) - std::pow(sigma_(1), 2)) +
  //           std::pow(sigma0, 2) * (1.0 / (1.0 - std::pow(eps, 2)));

  //   factor_int << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 / std::pow(eps, 2),
  //       qe(this->dir, 0) / qe(this->dir, 2);

  //   auto term_root = std::sqrt(qe(2, this->dir));
  //   auto term_rat = -((std::pow(pos_(1), 2)) / qe(this->dir, 0)) -
  //                   ((std::pow(pos_(2), 2)) / qe(this->dir, 1)) -
  //                   ((std::pow(pos_(3), 2)) / qe(this->dir, 2));
  //   return factor_int(this->dir, axis) * std::exp(term_rat) / term_root;
  // };

  // gsl_function_pp<decltype(fun)> Fp(fun);
  // gsl_function *F = static_cast<gsl_function *>(&Fp);
  // gsl_integration_workspace *work_ptr =
  // gsl_integration_workspace_alloc(1000); T abs_error = 1.0e-8; T rel_error
  // = 1.0e-8; T result; T error;

  // for (auto i = 0; i < 3; i++) {
  //   axis = i;
  //   gsl_integration_qag(F, kurt, 1.0, abs_error, rel_error, 1000,
  //                       GSL_INTEG_GAUSS31, work_ptr, &result, &error);
  //   E_(i + 1) = 2.0 * factor * (pos_(i + 1) / (sigma0 * cst::pi)) * result;
  // }
  // updateFields();
  // gsl_integration_workspace_free(work_ptr);
  quadv<T> E, B;
  return state_type2<T>{E, B};
}

template <class T>
state_type2<T> GaussianBunch<T>::transformFields(state_type2<T> fields) const {
  quadv<T> E{0.0, 0.0, 0.0, 0.0};
  quadv<T> B{0.0, 0.0, 0.0, 0.0};

  Eigen::Matrix<T, 4, 4> e_transform, b_transform;
  e_transform << 0.0, 0.0, 0.0, 0.0, 0.0, this->particle.getGamma(), 0.0, 0.0,
      0.0, 0.0, this->particle.getGamma(), 0.0, 0.0, 0.0, 0.0, 1.0;
  b_transform << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -(this->particle.getGamma() * this->particle.getBeta() / cst::sol), 0.0,
      0.0, this->particle.getGamma() * this->particle.getBeta() / cst::sol, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0;
  E = fields[0].transpose() * (e_transform);
  B = fields[0].transpose() * (b_transform);
  return state_type2<T>{E, B};
}

template <class T>
BunchEMField<T>::BunchEMField() : use_periodicity(true), local_time(.0) {}

template <class T>
void BunchEMField<T>::addBunch(std::unique_ptr<Bunch<T>> bunch) {
  bunches.push_back(std::move(bunch));
}

template <class T> void BunchEMField<T>::usePeriodicity(bool use) {
  use_periodicity = use;
}

template <class T> quadv<T> BunchEMField<T>::EfieldAt(quadv<T> quad) const {
  quadv<T> E{0.0, 0.0, 0.0, 0.0};
  for (auto &bunch : bunches) {
    if (use_periodicity) {
      auto tloc = (quad(0) / cst::sol);
      int rem = (int)std::floor(tloc / (bunch->getBunchPeriod()));
      tloc = tloc - rem * (1 * bunch->getBunchPeriod());
      if (tloc > (bunch->getBunchPeriod() / 2)) {
        quad(0) = (tloc - bunch->getBunchPeriod()) * cst::sol;
      } else {
        quad(0) = tloc * cst::sol;
      }
      E += bunch->EfieldAt(quad);
    } else {
      E += bunch->EfieldAt(quad);
    }
  }
  return E;
}

template <class T> quadv<T> BunchEMField<T>::MagfieldAt(quadv<T> quad) const {
  quadv<T> B;
  B << 0.0, 0.0, 0.0, 0.0;
  return B;
}

template <class T>
state_type2<T> BunchEMField<T>::EMfieldAt(quadv<T> quad) const {
  quadv<T> E{0.0, 0.0, 0.0, 0.0};
  quadv<T> B{0.0, 0.0, 0.0, 0.0};
  for (auto &bunch : bunches) {
    if (use_periodicity) {
      auto tloc = (quad(0) / cst::sol);
      int rem = (int)std::floor(tloc / (bunch->getBunchPeriod()));
      tloc = tloc - rem * (1 * bunch->getBunchPeriod());
      if (tloc > (bunch->getBunchPeriod() / 2)) {
        quad(0) = (tloc - bunch->getBunchPeriod()) * cst::sol;
      } else {
        quad(0) = tloc * cst::sol;
      }
      state_type2<T> temp_EM = bunch->EMfieldAt(quad);
      E += temp_EM[0];
      B += temp_EM[1];
    } else {
      state_type2<T> temp_EM = bunch->EMfieldAt(quad);
      E += temp_EM[0];
      B += temp_EM[1];
    }
  }
  return state_type2<T>{E, B};
} // namespace SpaceCharge

template class Bunch<double>;
template class GaussianBunch<double>;
template class BunchEMField<double>;

} // namespace SpaceCharge