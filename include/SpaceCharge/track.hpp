#ifndef TRACK_HPP
#define TRACK_HPP

#include "particle.hpp"
#include "point_cloud.hpp"

namespace SpaceCharge
{

  template <class T>
  class Track
  {
  private:
    /* data */
    Particle<T> particle;
    T t0;
    quadv<T> pos0;
    quadv<T> v0;

    std::vector<state_type> state;

  public:
    Track(Particle<T> particule, quadv);
    ~Track();
  };
} // namespace SpaceCharge

#endif