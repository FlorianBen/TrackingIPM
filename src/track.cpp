#include <boost/numeric/odeint.hpp>
#include <h5cpp/hdf5.hpp>
#include <thread>

#include "SpaceCharge/alogger.hpp"
#include "SpaceCharge/blosc_filter.h"
#include "SpaceCharge/definitions.hpp"
#include "SpaceCharge/track.hpp"
#include "SpaceCharge/track_h5.hpp"

namespace SpaceCharge {

template <class T>
Track<T>::Track()
    : particle("proton", 1, SpaceCharge::cst::mproton,
           SpaceCharge::cst::lfactor::beta, 0.5) {

           }

template <class T>
Track<T>::Track(Particle<T> part, quadv<T> pos0, quadv<T> v0,
                FieldSPS<T> &fieldmanager)
    : particle(part), pos0{pos0}, v0{v0}, fieldmanager(fieldmanager) {}

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
      pos.push_back(x[0]);
      // states.push_back(x);
    }
  };

  // Solve ODE
  integrate_const(stepper, lorentz, init, 0.0, 3e-9, 0.001e-9, observer);
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

template <class T>
void Track<T>::save(hdf5::node::Group group, const uint id) const {
  using namespace hdf5;
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  dcpl.chunk(hdf5::Dimensions{256});
  auto filter = std::make_unique<hdf5::filter::Deflate>(8u);
  filter->operator()(dcpl);

  node::Group root_group =
      group.create_group("track_" + std::to_string(id), lcpl);
  // node::Dataset dataset(root_group, "time",
  // datatype::create<std::vector<T>>(),
  //                       dataspace::create(times), lcpl, dcpl);

  // dataset.write(times);

  // auto type = datatype::create<SpaceCharge::quadv<T>>();
  auto dset = group.create_dataset(
      "pos", datatype::create<SpaceCharge::Track<double>>(),
      dataspace::create(*this), dcpl, lcpl);
  dset.write(*this);

  // VectorAppender<T> positions_appender(temp, "position");
  // for (auto state : states) {
  //   positions_appender(state[0]);
  // }
}

} // namespace SpaceCharge