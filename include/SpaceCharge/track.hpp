#ifndef TRACK_HPP
#define TRACK_HPP

#include "particle.hpp"

namespace SpaceCharge {

  template <class T> 
  class Track
  {
  private:
    /* data */
    Particle<T> particle;
    


  public:
    Track(/* args */);
    ~Track();
  };
  
  template <class T> 
  Track<T>::Track(/* args */)
  {
  }
  
  template <class T> 
  Track<T>::~Track()
  {
  }
  
}

#endif