#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "SpaceCharge/field/bunch.hpp"

namespace SpaceCharge {
using namespace dealii;

/**
 * @brief This class compute the field of a bunch by using FEM.
 *
 * @tparam dim
 */
template <int dim, class T> class FEMBunch : public Bunch<T> {
public:
  FEMBunch(Particle<T> particle, T ib, T dt, cst::dir dir = cst::dir::z);
  FEMBunch(int p_charge, T pmass, cst::lfactor factor, T lfactor_v, T ib, T dt,
           cst::dir dir = cst::dir::z);
  /**
   * @brief Compute the model.
   *
   */
  void run();

  /**
   * \brief Calculate the Electrical potential at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical potential (\f$V\f$)
   **/
  T potentialAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Electrical field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Electrical field vector (\f$V/m\f$)
   **/
  quadv<T> EfieldAt(quadv<T> quad) const override;
  /**
   * \brief Calculate the Magnetic field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return Magnetic field vector (\f$T\f$)
   **/
  quadv<T> MagfieldAt(quadv<T> quad) const override;

  /**
   * \brief Calculate the EM field at the given time and space.
   * \param[in] quad Time and space vector.
   * \return EM field double vectors.
   **/
  state_type2<T> EMfieldAt(quadv<T> quad) const override;

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;
  Triangulation<dim> triangulation;
  FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;
  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double> solution;
  Vector<double> system_rhs;
};
} // namespace SpaceCharge