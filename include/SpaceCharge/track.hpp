#ifndef TRACK_HPP
#define TRACK_HPP

#include <filesystem>

#include "particle.hpp"
#include "point_cloud.hpp"

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

public:
  Track(Particle<T> part, quadv<T> pos0, quadv<T> v0);

  void track();

  void save(std::filesystem::path file);
};

template class Track<double>;
template class Track<float>;

} // namespace SpaceCharge

#endif