#include "SpaceCharge/field/field_fem.hpp"
#include "SpaceCharge/core/alogger.hpp"
#include "SpaceCharge/core/definitions.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <Eigen/Eigen>

#include <fstream>
#include <iostream>

namespace SpaceCharge {

template <int dim> class RightHandSide : public Function<dim> {
protected:
  double current;
  double pulse_duration;
  double frequency;
  Eigen::Matrix<double, dim, 1> center;
  Eigen::Matrix<double, dim, dim> covar_m;

public:
  RightHandSide(Eigen::Matrix<double, dim, 1> center,
                Eigen::Matrix<double, dim, dim> covar)
      : center(center), covar_m(covar) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const {

  Eigen::Matrix<double, dim, 1> point_v;
  for (unsigned int i = 0; i < dim; ++i) {
    point_v(i) = p(i);
  }
  double Z_term = std::sqrt((2 * cst::pi * covar_m).determinant());
  double exp_term = std::exp(-0.5 * (point_v - center).transpose() *
                             covar_m.inverse() * (point_v - center));
  double f_value = (1.0 / Z_term) * exp_term;
  return f_value / cst::eps0;
}

template <int dim> class BoundaryValues : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double BoundaryValues<dim>::value(const Point<dim> &p,
                                  const unsigned int /*component*/) const {
  return 0;
}

template <int dim, class T>
FEMBunch<dim, T>::FEMBunch(Particle<T> particle, T ib, T dt, cst::dir dir)
    : Bunch<T>(particle, ib, dt, dir), fe(1), dof_handler(triangulation) {
  Logger::GetLogger()->info("FEMBunch: Create {}D FEMBunch", dim);
}

template <int dim, class T>
FEMBunch<dim, T>::FEMBunch(int p_charge, T pmass, cst::lfactor factor,
                           T lfactor_v, T ib, T dt, cst::dir dir)
    : Bunch<T>(p_charge, pmass, factor, lfactor_v, ib, dt, dir), fe(1),
      dof_handler(triangulation) {
  Logger::GetLogger()->info("FEMBunch: Create {}D FEMBunch", dim);
}

template <int dim, class T>
T FEMBunch<dim, T>::potentialAt(quadv<T> quad) const {
  T temp = 0.0;
  return temp;
}

template <int dim, class T>
quadv<T> FEMBunch<dim, T>::EfieldAt(quadv<T> quad) const {
  return EMfieldAt(quad)[0];
}

template <int dim, class T>
quadv<T> FEMBunch<dim, T>::MagfieldAt(quadv<T> quad) const {
  return EMfieldAt(quad)[1];
}

template <int dim, class T>
state_type2<T> FEMBunch<dim, T>::EMfieldAt(quadv<T> quad) const {
  quadv<T> E_{0.0, 0.0, 0.0, 0.0};
  quadv<T> B_{0.0, 0.0, 0.0, 0.0};
  return state_type2<T>{E_, B_};
}

template <int dim, class T> void FEMBunch<dim, T>::make_grid() {
  const Point<2> center(0, 0);
  const double inner_radius = 0.1, outer_radius = 0.125;
  GridGenerator::hyper_ball_balanced(triangulation, center, outer_radius);
  for (unsigned int step = 0; step < 7; ++step) {
    for (auto &cell : triangulation.active_cell_iterators()) {
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
        const double distance_from_center = center.distance(cell->vertex(v));
        if (distance_from_center < inner_radius) {
          cell->set_refine_flag();
          break;
        }
      }
    }
    triangulation.execute_coarsening_and_refinement();
  }
  Logger::GetLogger()->info("FEMBunch: {} total cells, {} active cells",
                            triangulation.n_active_cells(),
                            triangulation.n_cells());
}

// template <int dim, class T> void FEMBunch<3, T>::make_grid() {
//   GridGenerator::hyper_cube(triangulation, -1, 1);
//   triangulation.refine_global(4);
//   Logger::GetLogger()->info("FEMBunch: {} total cells, {} active cells",
//                             triangulation.n_active_cells(),
//                             triangulation.n_cells());
// }

template <int dim, class T> void FEMBunch<dim, T>::setup_system() {
  dof_handler.distribute_dofs(fe);
  Logger::GetLogger()->info("FEMBunch: {} degrees of freedom",
                            dof_handler.n_dofs());
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim, class T> void FEMBunch<dim, T>::assemble_system() {
  QGauss<dim> quadrature_formula(fe.degree + 1);
  Eigen::Matrix<double, dim, 1> center;
  center << .0, .0;
  Eigen::Matrix<double, dim, dim> covar;
  covar << .001, .0005, .0005, .001;

  RightHandSide<dim> right_hand_side(center, covar);
  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;
    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) +=
              (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
               fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
               fe_values.JxW(q_index));           // dx
        const auto &x_q = fe_values.quadrature_point(q_index);
        cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                        right_hand_side.value(x_q) *        // f(x_q)
                        fe_values.JxW(q_index));            // dx
      }
    cell->get_dof_indices(local_dof_indices);
    for (const unsigned int i : fe_values.dof_indices()) {
      for (const unsigned int j : fe_values.dof_indices())
        system_matrix.add(local_dof_indices[i], local_dof_indices[j],
                          cell_matrix(i, j));
      system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(
      dof_handler, 0, BoundaryValues<dim>(), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution,
                                     system_rhs);
}

template <int dim, class T> void FEMBunch<dim, T>::solve() {
  SolverControl solver_control(8000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  Logger::GetLogger()->info(
      "FEMBunch: {}  CG iterations needed to obtain convergence.",
      solver_control.last_step());
}

template <int dim, class T> void FEMBunch<dim, T>::output_results() const {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.add_data_vector(system_rhs, "rhs");
  data_out.build_patches();
  std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk");
  data_out.write_vtk(output);
}

template <int dim, class T> void FEMBunch<dim, T>::run() {
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}

template class FEMBunch<2, double>;
// template class FEMBunch<3, double>;

} // namespace SpaceCharge