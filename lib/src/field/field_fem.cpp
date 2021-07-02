#include "SpaceCharge/field/field_fem.hpp"

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

#include <fstream>
#include <iostream>

namespace SpaceCharge {

template <int dim> class RightHandSide : public Function<dim> {
protected:
  double current;
  double pulse_duration;
  double frequency;
  Point<dim> sigma;
  Point<dim> center;

public:
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;
};

template <int dim> class BoundaryValues : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;
};

template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const {
  double return_value = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 4.0 * std::pow(p(i), 4.0);
  return return_value;
}

template <int dim>
double BoundaryValues<dim>::value(const Point<dim> &p,
                                  const unsigned int /*component*/) const {
  return 0;
}

template <int dim>
FEMBunch<dim>::FEMBunch() : fe(1), dof_handler(triangulation) {}

template <> void FEMBunch<2>::make_grid() {
  const Point<2> center(0, 0);
  const double inner_radius = 0.050, outer_radius = 0.125;
  GridGenerator::hyper_ball_balanced(triangulation, center, outer_radius);
  for (unsigned int step = 0; step < 2; ++step) {
    for (auto &cell : triangulation.active_cell_iterators()) {
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v) {
        const double distance_from_center = center.distance(cell->vertex(v));
        // if (distance_from_center < inner_radius) {
        cell->set_refine_flag();
        // break;
        //}
      }
    }
    triangulation.execute_coarsening_and_refinement();
  }
}

template <> void FEMBunch<3>::make_grid() {
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(4);
  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}

template <int dim> void FEMBunch<dim>::setup_system() {
  dof_handler.distribute_dofs(fe);
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim> void FEMBunch<dim>::assemble_system() {
  QGauss<dim> quadrature_formula(fe.degree + 1);
  RightHandSide<dim> right_hand_side;
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

template <int dim> void FEMBunch<dim>::solve() {
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}

template <int dim> void FEMBunch<dim>::output_results() const {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.add_data_vector(system_rhs, "rhs");
  data_out.build_patches();
  std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk");
  data_out.write_vtk(output);
}

template <int dim> void FEMBunch<dim>::run() {
  std::cout << "Solving problem in " << dim << " space dimensions."
            << std::endl;
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}

template class FEMBunch<2>;
template class FEMBunch<3>;

} // namespace SpaceCharge