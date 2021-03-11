#ifndef TRACK_HPP
#define TRACK_HPP

#include <filesystem>
#include <h5cpp/hdf5.hpp>

#include "SpaceCharge/fields.hpp"
#include "SpaceCharge/particle.hpp"
#include "SpaceCharge/point_cloud.hpp"

namespace SpaceCharge {

template <class T> class Track {
private:
  /* data */
  Particle<T> particle;
  T t0;
  quadv<T> pos0;
  quadv<T> v0;

  std::vector<state_type2<T>> states;
  std::vector<quadv<T>> pos;
  std::vector<T> times;

  FieldSPS<T> fieldmanager;

  bool filter(const quadv<T> &pos) const;

public:
  Track();
  Track(Particle<T> part, quadv<T> pos0, quadv<T> v0,
        FieldSPS<T> &fieldmanager);

  void track();

  quadv<T> *data();

  const quadv<T> *data() const;

  size_t size() const;

  void save(hdf5::node::Group group, const uint id) const;
};

template class Track<double>;
template class Track<float>;

} // namespace SpaceCharge

#endif