#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>

#include "SpaceCharge/field/field_map.hpp"

namespace SpaceCharge {

// TODO: Refactor to Grid.
// TODO: Split into different class like GridCreator etc...

template <typename T>
FieldMap<T>::FieldMap()
    : qsize(quadv<size_t>{0, 0, 0, 0}), qstep(quadv<T>{.0, .0, .0, .0}),
      qoffset(quadv<T>{.0, .0, .0, .0}), time(.0), data_(0) {}

template <typename T>
FieldMap<T>::FieldMap(quadv<size_t> size, quadv<T> step, quadv<T> offset,
                      T time)
    : qsize(size), qstep(step), qoffset(offset), time(.0),
      data_(qsize(1) * qsize(2) * qsize(3)) {}

template <typename T> size_t FieldMap<T>::size() const { return data_.size(); }

template <typename T> quadv<size_t> FieldMap<T>::sizes() const { return qsize; }

template <typename T> quadv<T> FieldMap<T>::steps() const { return qstep; }

template <typename T> quadv<T> FieldMap<T>::offsets() const { return qoffset; }

template <typename T> size_t FieldMap<T>::nx() const { return qsize(1); }

template <typename T> size_t FieldMap<T>::ny() const { return qsize(2); }

template <typename T> size_t FieldMap<T>::nz() const { return qsize(3); }

template <typename T> quadv<T> *FieldMap<T>::data() { return data_.data(); }

template <typename T> const quadv<T> *FieldMap<T>::data() const {
  return data_.data();
}

template <typename T>
const std::vector<quadv<T>> &FieldMap<T>::getVector() const {
  return data_;
}

template <typename T> std::vector<quadv<T>> &FieldMap<T>::getVector() {
  return data_;
}

template <typename T>
quadv<T> FieldMap<T>::operator()(size_t x, size_t y, size_t z) const {
  return data_[z + x * nz() + y * nx() * nz()];
}

template <typename T>
quadv<T> &FieldMap<T>::operator()(size_t x, size_t y, size_t z) {
  return data_[z + x * nz() + y * nx() * nz()];
}

template <typename T> void FieldMapInterpolate<T>::addField(FieldSP<T> &field) {
  input_fields.push_back(std::move(field));
}

template <typename T> void FieldMapInterpolate<T>::computeField() {
  tbb::parallel_for(
      tbb::blocked_range3d<int>(0, FieldMap<T>::nz(), 0, FieldMap<T>::ny(), 0,
                                FieldMap<T>::nx()),
      [&](const tbb::blocked_range3d<int> &r) {
        for (int y = r.rows().begin(), y_end = r.rows().end(); y < y_end; y++) {
          for (int x = r.cols().begin(), x_end = r.cols().end(); x < x_end;
               x++) {
            for (int z = r.pages().begin(), z_end = r.pages().end(); z < z_end;
                 z++) {
              for (auto &field : input_fields) {
                quadv<T> pos{
                    FieldMap<T>::qoffset(0),
                    x * FieldMap<T>::qstep(1) + FieldMap<T>::qoffset(1),
                    y * FieldMap<T>::qstep(2) + FieldMap<T>::qoffset(2),
                    z * FieldMap<T>::qstep(3) + FieldMap<T>::qoffset(3)};
                auto res = field->EMfieldAt(pos);
                FieldMap<T>::data_[z + x * FieldMap<T>::nz() +
                                   y * FieldMap<T>::nx() * FieldMap<T>::nz()] =
                    res[0];
                // FieldMap<T>::data_B[z + x * FieldMap<T>::nz() +
                //                     y * FieldMap<T>::nx() * FieldMap<T>::nz()] =
                //     res[1];
              }
            }
          }
        }
      });
}

template class FieldMap<double>;
template class FieldMapInterpolate<double>;

} // namespace SpaceCharge