#include "SpaceCharge/fields.hpp"

#include <iostream>

#include "SpaceCharge/csv.hpp"

namespace SpaceCharge {
template <class T> Field<T>::Field() {}

template <class T> Field<T>::~Field() {}

template <class T>
FieldBunch<T>::FieldBunch() : use_periodicity(true), local_time(.0) {}

template <class T>
void FieldBunch<T>::addBunch(std::unique_ptr<Bunch<T>> bunch) {
  bunches.push_back(std::move(bunch));
}

template <class T> void FieldBunch<T>::usePeriodicity(bool use) {
  use_periodicity = use;
}

template <class T> Eigen::Matrix<T, 4, 1> FieldBunch<T>::EfieldAt(quadv quad) {
  quadv E;
  quadv pos1;
  quadv pos2;
  quadv pos3;
  E << 0.0, 0.0, 0.0, 0.0;
  for (auto &bunch : bunches) {
    if (use_periodicity) {
      auto tloc = (quad(0) / cst::sol);
      int rem = (int)std::floor(tloc / (bunch->getBunchPeriod()));
      tloc = tloc - rem * (1 * bunch->getBunchPeriod());
      if (tloc > (bunch->getBunchPeriod() / 2)) {
        quad(0) = (tloc - bunch->getBunchPeriod()) * cst::sol;
      } else {
        quad(0) = tloc * cst::sol;
      }

      E += bunch->EfieldAt(quad);
    } else {
      E += bunch->EfieldAt(quad);
    }
  }

  return E;
}

template <class T>
Eigen::Matrix<T, 4, 1> FieldBunch<T>::MagfieldAt(quadv quad) {

  return quad;
}

template <class T> FieldCOMSOL<T>::FieldCOMSOL() {}

template <class T> FieldCOMSOL<T>::~FieldCOMSOL() {}

template <class T> Eigen::Matrix<T, 4, 1> FieldCOMSOL<T>::EfieldAt(quadv quad) {
  return quad;
}

template <class T>
Eigen::Matrix<T, 4, 1> FieldCOMSOL<T>::MagfieldAt(quadv quad) {
  return quad;
}

template <class T>
void FieldCOMSOL<T>::loadEfield(const std::string filename, const quadv offset,
                                const double scale) {
  constexpr auto ncols = 6;
  io::CSVReader<ncols> in(filename);
  quadv tmp_pos;
  tmp_pos << .0, .0, .0, .0;
  quadv tmp_val;
  tmp_val << .0, .0, .0, .0;
  while (in.read_row(tmp_pos(1), tmp_pos(2), tmp_pos(3), tmp_val(1), tmp_val(2),
                     tmp_val(3))) {
    tmp_pos -= offset;
    pointcloud_posE.pts.push_back(tmp_pos);
    fieldE.pts.push_back(tmp_val * scale);
  }
}

template <class T> void FieldCOMSOL<T>::create_Eindex(const int leaf_size) {
  index = new kd_tree_nanoflann(
      3, pointcloud_posE, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size));
  index->buildIndex();
}

template <class T>
size_t FieldCOMSOL<T>::interpolateNN(const quadv pos, quadv &fieldv,
                                     size_t size) const {
  T query_pt[4] = {pos[0], pos[1], pos[2], pos[3]};
  size_t nb_result = size;
  std::vector<size_t> ret_index(nb_result);
  std::vector<T> out_dist_sqr(nb_result);
  nb_result = index->knnSearch(&query_pt[0], nb_result, &ret_index[0],
                               &out_dist_sqr[0]);

  fieldv = fieldE.pts[ret_index[0]];
  return (nb_result > 0) ? nb_result : -1;
}

template <class T>
size_t FieldCOMSOL<T>::interpolateNN2(const quadv pos, quadv &fieldv,
                                      size_t size) const {
  T query_pt[4] = {pos[0], pos[1], pos[2], pos[3]};
  const size_t nb_results = 1;
  size_t ret_index;
  T out_dist_sqr;
  nanoflann::KNNResultSet<T> resultSet(nb_results);
  resultSet.init(&ret_index, &out_dist_sqr);
  index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(size));
  fieldv = fieldE.pts[ret_index];
  return (nb_results > 0) ? nb_results : -1;
}

template <class T>
void FieldCOMSOL<T>::interpolateRBF(
    const quadv pos, quadv &fieldv, const int nbNN, const float order,
    const std::function<T(const T, const T)> kernel) const {
  T query_pt[4] = {pos[0], pos[1], pos[2], pos[3]};

  std::vector<size_t> ret_index(nbNN);
  std::vector<T> out_dist_sqr(nbNN);
  auto nb_result =
      index->knnSearch(&query_pt[0], nbNN, &ret_index[0], &out_dist_sqr[0]);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matW(nb_result, nb_result);

  for (auto i = 0; i < nb_result; i++) {
    for (auto j = 0; j < nb_result; j++) {
      matW(i, j) = kernel(distanceFrom(pointcloud_posE.pts[ret_index[i]],
                                       pointcloud_posE.pts[ret_index[j]]),
                          order);
    }
  }
  Eigen::Matrix<T, Eigen::Dynamic, 1> Wkx(nb_result), Wky(nb_result),
      Wkz(nb_result), Fkx(nb_result), Fky(nb_result), Fkz(nb_result);
  for (auto ind = 0; ind < ret_index.size(); ind++) {
    Fkx(ind) = fieldE.pts[ret_index[ind]](1);
    Fky(ind) = fieldE.pts[ret_index[ind]](2);
    Fkz(ind) = fieldE.pts[ret_index[ind]](3);
  }
  Wkx = matW.colPivHouseholderQr().solve(Fkx);
  Wky = matW.colPivHouseholderQr().solve(Fky);
  Wkz = matW.colPivHouseholderQr().solve(Fkz);
  T Fix = .0, Fiy = .0, Fiz = .0;
  for (auto ind = 0; ind < ret_index.size(); ind++) {
    auto distance = distanceFrom(pos, pointcloud_posE.pts[ret_index[ind]]);
    Fix = Fix + kernel(distance, order) * Wkx(ind);
    Fiy = Fiy + kernel(distance, order) * Wky(ind);
    Fiz = Fiz + kernel(distance, order) * Wkz(ind);
  }
  fieldv(1) = Fix;
  fieldv(2) = Fiy;
  fieldv(3) = Fiz;
}

template class Field<double>;
template class FieldBunch<double>;
template class FieldCOMSOL<double>;
template class FieldCOMSOL<float>;

}; // namespace SpaceCharge