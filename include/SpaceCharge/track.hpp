#ifndef TRACK_HPP
#define TRACK_HPP

#include <filesystem>

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
  std::vector<T> times;

  FieldSPS<T> fieldmanager;

public:
  Track(Particle<T> part, quadv<T> pos0, quadv<T> v0, FieldSPS<T> &fieldmanager);

  void track();

  void save(std::filesystem::path file);
};

template class Track<double>;
template class Track<float>;

} // namespace SpaceCharge

#endif