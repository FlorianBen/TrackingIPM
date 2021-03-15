#include <iostream>

#include "SpaceCharge/alogger.hpp"
#include "SpaceCharge/csv.hpp"
#include "SpaceCharge/fields.hpp"

// TODO: Implémenter les fonctions

namespace SpaceCharge {
template <class T> EMField<T>::EMField() {}

template <class T> EMField<T>::~EMField() {}

template <class T>
ConstantEMField<T>::ConstantEMField(state_type2<T> emfield)
    : emfield(emfield) {}

template <class T> quadv<T> ConstantEMField<T>::EfieldAt(quadv<T> quad) const {
  return emfield[0];
}

template <class T>
quadv<T> ConstantEMField<T>::MagfieldAt(quadv<T> quad) const {
  return emfield[1];
}

template <class T>
state_type2<T> ConstantEMField<T>::EMfieldAt(quadv<T> quad) const {
  return emfield;
}

template <class T> CSVFileEMField<T>::CSVFileEMField() {}

template <class T> CSVFileEMField<T>::~CSVFileEMField() {}

template <class T> quadv<T> CSVFileEMField<T>::EfieldAt(quadv<T> quad) const {
  quadv<T> E;
  E << 0.0, 0.0, 0.0, 0.0;
  interpolateRBF(quad, E, 7, 1.0, kernel_exp);
  return E;
}

template <class T> quadv<T> CSVFileEMField<T>::MagfieldAt(quadv<T> quad) const {
  quadv<T> B;
  B << 0.0, 0.0, 0.0, 0.0;
  return B;
}

template <class T>
state_type2<T> CSVFileEMField<T>::EMfieldAt(quadv<T> quad) const {
  return state_type2<T>{EfieldAt(quad), MagfieldAt(quad)};
}

template <class T>
void CSVFileEMField<T>::loadEfield(const std::string filename,
                                   const quadv<T> offset, const double scale,
                                   const int leaf_size) {
  SC_INFO("FieldCOMSOL: Load file {}", filename);
  constexpr auto ncols = 6;
  io::CSVReader<ncols> in(filename);
  quadv<T> tmp_pos;
  tmp_pos << .0, .0, .0, .0;
  quadv<T> tmp_val;
  tmp_val << .0, .0, .0, .0;
  while (in.read_row(tmp_pos(1), tmp_pos(2), tmp_pos(3), tmp_val(1), tmp_val(2),
                     tmp_val(3))) {
    tmp_pos -= offset;
    pointcloud_posE.pts.push_back(tmp_pos);
    fieldE.pts.push_back(tmp_val * scale);
  }
  SC_INFO("FieldCOMSOL: File loaded");
  create_Eindex(leaf_size);
}

template <class T> void CSVFileEMField<T>::create_Eindex(const int leaf_size) {
  index = new kd_tree_nanoflann(
      3, pointcloud_posE, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size));
  index->buildIndex();
  SC_INFO("FieldCOMSOL: Index created");
}

template <class T>
size_t CSVFileEMField<T>::interpolateNN(const quadv<T> pos, quadv<T> &fieldv,
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
size_t CSVFileEMField<T>::interpolateNN2(const quadv<T> pos, quadv<T> &fieldv,
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
void CSVFileEMField<T>::interpolateRBF(
    const quadv<T> pos, quadv<T> &fieldv, const int nbNN, const float order,
    const std::function<T(const T, const T)> kernel) const {
  T query_pt[4] = {pos[0], pos[1], pos[2], pos[3]};

  std::vector<size_t> ret_index(nbNN);
  std::vector<T> out_dist_sqr(nbNN);
  auto nb_result =
      index->knnSearch(&query_pt[0], nbNN, &ret_index[0], &out_dist_sqr[0]);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matW(nb_result, nb_result);

  for (auto i = 0; i < nb_result; i++) {
    for (auto j = 0; j < nb_result; j++) {
      matW(i, j) = kernel(distanceFrom<T>(pointcloud_posE.pts[ret_index[i]],
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
    auto distance = distanceFrom<T>(pos, pointcloud_posE.pts[ret_index[ind]]);
    Fix = Fix + kernel(distance, order) * Wkx(ind);
    Fiy = Fiy + kernel(distance, order) * Wky(ind);
    Fiz = Fiz + kernel(distance, order) * Wkz(ind);
  }
  fieldv(1) = Fix;
  fieldv(2) = Fiy;
  fieldv(3) = Fiz;
}

template <class T> EMFieldsManager<T>::EMFieldsManager() {
  // TODO: Implement
}

template <class T> void EMFieldsManager<T>::addField(FieldSP<T> &field) {
  fields.push_back(std::move(field));
}

template <class T> quadv<T> EMFieldsManager<T>::EfieldAt(quadv<T> quad) const {
  quadv<T> E;
  E << 0.0, 0.0, 0.0, 0.0;
  for (auto &field : fields) {
    E += field->EfieldAt(quad);
  }
  return E;
}

template <class T>
quadv<T> EMFieldsManager<T>::MagfieldAt(quadv<T> quad) const {
  quadv<T> B;
  B << 0.0, 0.0, 0.0, 0.0;
  for (auto &field : fields) {
    B += field->MagfieldAt(quad);
  }
  return B;
}

template <class T>
state_type2<T> EMFieldsManager<T>::EMfieldAt(quadv<T> quad) const {
  quadv<T> E;
  quadv<T> B;
  E << 0.0, 0.0, 0.0, 0.0;
  B << 0.0, 0.0, 0.0, 0.0;
  for (auto &field : fields) {
    auto em = field->EMfieldAt(quad);
    E += em[0];
    B += em[1];
  }
  return state_type2<T>{E, B};
}

template class EMField<double>;
template class EMField<float>;
template class ConstantEMField<double>;
template class ConstantEMField<float>;
template class CSVFileEMField<double>;
template class CSVFileEMField<float>;
template class EMFieldsManager<double>;
template class EMFieldsManager<float>;

}; // namespace SpaceCharge