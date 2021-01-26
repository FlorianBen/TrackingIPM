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
Track<T>::Track(Particle<T> part, quadv<T> pos0, quadv<T> v0,
                FieldSPS<T> &fieldmanager)
    : particle(part), pos0{pos0}, v0{v0}, fieldmanager(fieldmanager) {}

template <class T> void Track<T>::track() {
  using namespace boost::numeric::odeint;
  state_type2<T> init{pos0, v0};
  runge_kutta4<state_type2<T>> stepper;

  // Define lorentz equation
  auto lorentz = [&](const state_type2<T> &x, state_type2<T> &dxdt,
                     const double t) {
    quadv<T> Efield, Bfield;
    state_type2<T> EMfield = this->fieldmanager->EMfielddAt(x[0]);
    dxdt[0] = x[1];
    dxdt[1] =
        ((particle.getCharge() * SpaceCharge::cst::e) /
         (particle.getMass() * sqrt(1 - SpaceCharge::scalar_prod(x[1], x[1]) /
                                            (cst::sol * cst::sol)))) *
        (Efield + SpaceCharge::vect_prod(x[1], Bfield));
  };

  // Define observer
  auto observer = [&](state_type2<T> &x, T t) {
    times.push_back(t);
    x[0](0) = t;
    states.push_back(x);
  };

  // Solve ODE
  integrate_const(stepper, lorentz, init, 0.0, 3e-9, 0.001e-9, observer);
}

template <class T>
void Track<T>::save(hdf5::node::Group group, const uint id) const {
  using namespace hdf5;
  hdf5::property::LinkCreationList lcpl;
  hdf5::property::DatasetCreationList dcpl;
  dcpl.layout(hdf5::property::DatasetLayout::CHUNKED);
  dcpl.chunk(hdf5::Dimensions{256});
  //auto filter = std::make_unique<hdf5::filter::Deflate>(8u);
  //filter->operator()(dcpl);

  node::Group root_group =
      group.create_group("track_" + std::to_string(id), lcpl);
  node::Dataset dataset(root_group, "time", datatype::create<std::vector<T>>(),
                        dataspace::create(times), lcpl, dcpl);

  dataset.write(times);

  auto type = datatype::create<SpaceCharge::quadv<T>>();

  hdf5::dataspace::Simple dataspace({0}, {hdf5::dataspace::Simple::UNLIMITED});
  auto temp =
      hdf5::node::Dataset(root_group, "positions", type, dataspace, lcpl, dcpl);
  VectorAppender<T> positions_appender(temp, "position");
  for (auto state : states) {
    positions_appender(state[0]);
  }
}

} // namespace SpaceCharge