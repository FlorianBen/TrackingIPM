#include "SpaceCharge/bunch.hpp"

namespace SpaceCharge {

template <class T> Bunch<T>::Bunch(){};

template <class T>
Bunch<T>::Bunch(std::string pdg_name, T ib, T dt)
    : ib(ib), dt(dt){
        check_particle_def(pdg_name);
              };

//template class Bunch<double>;
//template class Bunch<float>;

} // namespace SpaceCharge