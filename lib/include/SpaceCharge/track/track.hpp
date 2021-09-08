#ifndef TRACK_HPP
#define TRACK_HPP

#include <filesystem>

#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/core/point_cloud.hpp"
#include "SpaceCharge/field/fields.hpp"

namespace SpaceCharge {

/**
 * \class Track track.hpp
 * \brief Track class describes a particule trajectory in a field.
 * The particule and field can be defined, then the trjaectory is computed by
 * mean of ODE.
 **/
template <class T> class Track {
private:
  /* data */
  Particle<T> particle;
  T t0;
  quadv<T> pos0;
  quadv<T> v0;

  std::vector<quadv<T>> pos, speed;
  std::vector<T> times;

  FieldSPS<T> fieldmanager;

  bool filter(const quadv<T> &pos) const;

public:
  /**
   * \brief Default constructor.
   * Construct an empty Track.
   **/
  Track();

  /**
   * \brief Constructor.
   * \param[in] part E and B vectors.
   * \param[in] pos0 E and B vectors.
   * \param[in] v0 E and B vectors.
   * \param[in] fieldmanager E and B vectors.
   **/
  Track(Particle<T> part, quadv<T> pos0, quadv<T> v0,
        FieldSPS<T> &fieldmanager);

  /**
   * \brief Run the ODE tracking.
   **/
  void track();

  /**
   * \brief Run the boris tracking.
   **/
  void track2();

  /**
   * \brief Get the raw data vector.
   * \return data pointer.
   **/
  quadv<T> *data();

  /**
   * \brief Get the const raw data vector.
   * \return data pointer.
   **/
  const quadv<T> *data() const;

  const std::vector<quadv<T>> &getPosVector() const;

  const std::vector<quadv<T>> &getSpeedVector() const;

  /**
   * \brief Get the size of the track.
   * \return Track size.
   **/
  size_t size() const;
};

template class Track<double>;
template class Track<float>;

} // namespace SpaceCharge

#endif