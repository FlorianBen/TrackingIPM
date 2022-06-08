#include <boost/numeric/odeint.hpp>
#include <thread>

#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/definitions.hpp"
#include "SpaceCharge/io/blosc_filter.h"
#include "SpaceCharge/track/track.hpp"

namespace SpaceCharge {

// TODO: Add custom stepper.
// TODO: Split into a generic Track (without ODE)
// TODO: Add stepper options or configuration.

template <class T>
Track<T>::Track()
    : particle("proton", 1, SpaceCharge::cst::mproton,
               SpaceCharge::cst::lfactor::beta, 0.5) {}

template <class T>
Track<T>::Track(Particle<T> part, quadv<T> pos0, quadv<T> v0,
                FieldSPS<T> &fieldmanager)
    : particle(part), pos0{pos0}, v0{v0}, fieldmanager(fieldmanager) {}

template <class T> Track<T>::~Track() {}


template <class T> void Track<T>::track() {
  using namespace boost::numeric::odeint;
  state_type2<T> init{pos0, v0};
  runge_kutta4<state_type2<T>> stepper;

  // Define lorentz equation
  auto lorentz = [&](state_type2<T> &x, state_type2<T> &dxdt, const double t) {
    x[0](0) = SpaceCharge::cst::sol * t;
    state_type2<T> EMfield = this->fieldmanager->EMfieldAt(x[0]);
    dxdt[0] = x[1];
    dxdt[1] =
        ((particle.getCharge() * SpaceCharge::cst::e) /
         (particle.getMass() * sqrt(1 - SpaceCharge::scalar_prod(x[1], x[1]) /
                                            (cst::sol * cst::sol)))) *
        (EMfield[0] + SpaceCharge::vect_prod(x[1], EMfield[1]));
  };

  // Define observer
  auto observer = [&](state_type2<T> &x, T t) {
    if (filter(x[0])) {
      times.push_back(t);
      x[0](0) = t;
      x[1](0) = t;
      pos.push_back(x[0]);
      speed.push_back(x[1]);
    }
  };

  // Solve ODE
  integrate_const(stepper, lorentz, init, 0.0, 100e-9, 0.01e-9, observer);

  auto reduce = true;
  if(reduce) {
    auto posi_red = *std::begin(pos);
    auto posf_red = *std::rbegin(pos);
    auto speedi_red = *std::begin(speed);
    auto speedf_red = *std::rbegin(speed);
    pos.clear();
    pos.push_back(posi_red);
    pos.push_back(posf_red);
    speed.clear();
    speed.push_back(speedi_red);
    speed.push_back(speedf_red);
    
  }
}

template <class T> void Track<T>::track2(){};

template <typename T>
const std::vector<quadv<T>> &Track<T>::getPosVector() const {
  return pos;
}

template <typename T>
const std::vector<quadv<T>> &Track<T>::getSpeedVector() const {
  return speed;
}

template <typename T> quadv<T> *Track<T>::data() { return pos.data(); }

template <typename T> const quadv<T> *Track<T>::data() const {
  return pos.data();
}

template <typename T> size_t Track<T>::size() const { return pos.size(); }

template <class T> bool Track<T>::filter(const quadv<T> &pos) const {
  if ((pos(1) > 0.05) || (pos(1) < -0.05))
    return false;
  if ((pos(2) > 0.05) || (pos(2) < -0.05))
    return false;
  if ((pos(3) > 0.05) || (pos(3) < -0.05))
    return false;
  return true;
}

template class Track<double>;
template class Track<float>;

} // namespace SpaceCharge