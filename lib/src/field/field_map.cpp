#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>

#include "SpaceCharge/field/field_map.hpp"

namespace SpaceCharge {

template <typename T>
FieldMap<T>::FieldMap()
    : qsize(quadv<size_t>{0, 0, 0, 0}), qstep(quadv<T>{.0, .0, .0, .0}),
      qoffset(quadv<T>{.0, .0, .0, .0}), time(.0), data_(0) {}

template <typename T>
FieldMap<T>::FieldMap(quadv<size_t> size, quadv<T> step, quadv<T> offset, T time)
    : qsize(size), qstep(step), qoffset(offset),
      data_(qsize(1) * qsize(2) * qsize(3)), time(.0) {

      };

template <typename T> size_t FieldMap<T>::size() const { return data_.size(); }

template <typename T> size_t FieldMap<T>::nx() const { return qsize(1); }

template <typename T> size_t FieldMap<T>::ny() const { return qsize(2); }

template <typename T> size_t FieldMap<T>::nz() const { return qsize(3); }

template <typename T> quadv<T> *FieldMap<T>::data() { return data_.data(); }

template <typename T> const quadv<T> *FieldMap<T>::data() const {
  return data_.data();
}

template <typename T>
quadv<T> FieldMap<T>::operator()(size_t x, size_t y, size_t z) const {
  return data_[z + x * nz() + y * nx() * nz()];
}

template <typename T>
quadv<T> &FieldMap<T>::operator()(size_t x, size_t y, size_t z) {
  return data_[z + x * nz() + y * nx() * nz()];
}

template <typename T> void FieldMap<T>::addField(FieldSP<T> &field) {
  input_fields.push_back(std::move(field));
}

template <typename T> void FieldMap<T>::computeField() {
  tbb::parallel_for(
      tbb::blocked_range3d<int>(0, nz(), 0, ny(), 0, nx()),
      [&](const tbb::blocked_range3d<int> &r) {
        for (int y = r.rows().begin(), y_end = r.rows().end(); y < y_end; y++) {
          for (int x = r.cols().begin(), x_end = r.cols().end(); x < x_end;
               x++) {
            for (int z = r.pages().begin(), z_end = r.pages().end(); z < z_end;
                 z++) {
              for (auto &field : input_fields) {
                quadv<T> pos{qoffset(0), x * qstep(1) + qoffset(1),
                             y * qstep(2) + qoffset(2),
                             z * qstep(3) + qoffset(3)};
                data_[z + x * nz() + y * nx() * nz()] =
                    field->EMfieldAt(pos)[0];
              }
            }
          }
        }
      });
}

template class FieldMap<double>;

} // namespace SpaceCharge