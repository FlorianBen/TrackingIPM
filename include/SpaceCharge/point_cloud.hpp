#ifndef POINT_CLOUD_HPP
#define POINT_CLOUD_HPP

#undef Success // Conflit entre les libs 
#include <Eigen/Dense>

#include <iostream>

namespace SpaceCharge {

template <typename T> using quadv = Eigen::Matrix<T, 4, 1>;

template <typename T> using state_type2 = std::vector<Eigen::Matrix<T, 4, 1>>;

typedef std::vector<Eigen::Matrix<double, 4, 1>> state_type;

template <class T>
T distanceFrom(const Eigen::Matrix<T, 4, 1> &p1,
               const Eigen::Matrix<T, 4, 1> &p2) {
  T distsqr = 0;
  for (size_t i = 1; i < 4; ++i) {
    distsqr = distsqr + pow((p1(i) - p2(i)), 2);
  }
  return sqrt(distsqr);
}

/**
 * @brief Scalar/Dot product between two vector (3D).
 * @param[in] p1
 * @param[in] p2
 * @return The result of dot product.
 */
template <class T>
T scalar_prod(const Eigen::Matrix<T, 4, 1> &p1,
              const Eigen::Matrix<T, 4, 1> &p2) {
  T tmp = 0.0;
  for (size_t i = 1; i < 4; ++i)
    tmp += p1(i) * p2(i);
  return tmp;
}

/**
 * @brief Vector/Cross product between two vector (3D).
 * @param[in] p1
 * @param[in] p2
 * @return The result of cross product.
 */
template <class T>
Eigen::Matrix<T, 4, 1> vect_prod(const Eigen::Matrix<T, 4, 1> &p1,
                                 const Eigen::Matrix<T, 4, 1> &p2) {
  Eigen::Matrix<T, 4, 1> tmp;
  tmp << .0, .0, .0, .0;
  tmp(1) = p1(2) * p2(3) - p1(3) * p2(2);
  tmp(2) = p1(3) * p2(1) - p1(1) * p2(3);
  tmp(3) = p1(1) * p2(2) - p1(2) * p2(1);
  return tmp;
}

template <class T> struct PointCloud {
  /**
   * @struct point point_type.hpp "include/point_type.hpp"
   * @brief Adaptator class between nanoflann and point_type.
   */

  std::vector<quadv<T>> pts;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  // Returns the distance between the vector "p1[0:size-1]" and the data point
  // with index "idx_p2" stored in the class:
  inline T kdtree_distance(const T *p1, const size_t idx_p2,
                           size_t /*size*/) const {
    // const T d0 = p1(0) - pts[idx_p2](0);
    const T d1 = p1(1) - pts[idx_p2](1);
    const T d2 = p1(2) - pts[idx_p2](2);
    const T d3 = p1(3) - pts[idx_p2](3);
    return d1 * d1 + d2 * d2 + d3 * d3; // d0 * d0 +
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline T kdtree_get_pt(const size_t idx, int dim) const {
    if (dim == 0)
      return pts[idx](0);
    else if (dim == 1)
      return pts[idx](1);
    else if (dim == 2)
      return pts[idx](2);
    else
      return pts[idx](3);
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3
  //   for point clouds)
  template <class BBOX> bool kdtree_get_bbox(BBOX & /* bb */) const {
    return false;
  }
};

}; // namespace SpaceCharge

#endif