/*
 * ThermalProblem.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef THERMALPROBLEM_H_
#define THERMALPROBLEM_H_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/iterative_inverse.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>

#include "RightHandSide.h"
#include "BoundaryValues.h"
#include "BoundaryCondition.h"
#include "TemperatureBoundary.h"

using namespace dealii;

template <int dim>
class ThermalProblem
{

public:
	ThermalProblem () : fe(1), dof_handler(triangulation) {};
	~ThermalProblem() {};

	void run ();

private:
	void make_grid ();
	void setup_system();
	void assemble_system ();
	void solve ();
	void output_results () const;

	Triangulation<dim>   triangulation;
	FE_Q<dim>            fe;
	DoFHandler<dim>      dof_handler;
	BoundaryValues<dim>	 boundary_values;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double>       solution;
	Vector<double>       system_rhs;

};

// Public method: run
template<int dim>
void ThermalProblem<dim>::run()
{
	std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;

	make_grid();
	setup_system ();
	assemble_system ();
	solve ();
	output_results ();
}

// Private method: make_grid
template<int dim>
void ThermalProblem<dim>::make_grid()
{
	// For now just generate cube
	// Later include functionality to read in a mesh file?
	GridGenerator::hyper_cube(triangulation, -1, 1);
	triangulation.refine_global(4);

	// Create boundary conditions
	Point<dim> p1(true), p2(true);
	if (dim == 2) {
		p1(0) = -1;
		p1(1) = -1;
		p2(0) = 1;
		p2(1) = 1;
	}
	else if (dim ==3) {
		p1(0) = -1;
		p1(1) = -1;
		p1(2) = -1;
		p2(0) = 1;
		p2(1) = 1;
		p2(2) = 1;
	}
	else {
		Assert(false, ExcNotImplemented())
	}
	// For all boundaries in the rectangular bounds created by p1 and p2,
	// apply a specified temperature of 100
	BoundaryCondition<dim> * bc_p = new TemperatureBoundary<dim>(p1, p2, 100);
	boundary_values.add_boundary_condition(bc_p);
}

// Private method: setup_system
template<int dim>
void ThermalProblem<dim>::setup_system()
{
	dof_handler.distribute_dofs (fe);

	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << "\n";

	CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);

	system_matrix.reinit (sparsity_pattern);

	solution.reinit (dof_handler.n_dofs());
	system_rhs.reinit (dof_handler.n_dofs());
}

// Private method: assemble_system
template<int dim>
void ThermalProblem<dim>::assemble_system()
{

	QGauss<dim>  quadrature_formula(2);

	const RightHandSide<dim> right_hand_side;

	FEValues<dim> fe_values (fe, quadrature_formula,
	update_values   | update_gradients |
	update_quadrature_points | update_JxW_values);

	const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	const unsigned int   n_q_points    = quadrature_formula.size();

	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	for (; cell!=endc; ++cell)
	{
		fe_values.reinit (cell);
		cell_matrix = 0;
		cell_rhs = 0;

		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		{
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
					cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
							fe_values.shape_grad (j, q_point) *
							fe_values.JxW (q_point));

				cell_rhs(i) += (fe_values.shape_value (i, q_point) *
							right_hand_side.value (fe_values.quadrature_point (q_point)) *
							fe_values.JxW (q_point));
			}
		}

		cell->get_dof_indices (local_dof_indices);
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			for (unsigned int j=0; j<dofs_per_cell; ++j)
				system_matrix.add (local_dof_indices[i],
							local_dof_indices[j],
							cell_matrix(i,j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}


	std::map<unsigned int,double> boundary_values_map;
	VectorTools::interpolate_boundary_values (dof_handler,
				0,
				boundary_values,
				boundary_values_map);

	MatrixTools::apply_boundary_values (boundary_values_map,
		  system_matrix,
		  solution,
		  system_rhs);

}

// Private method: solve
template<int dim>
void ThermalProblem<dim>::solve()
{
	SolverControl           solver_control (1000, 1e-12);
	SolverCG<>              solver (solver_control);
	solver.solve (system_matrix, solution, system_rhs,
				PreconditionIdentity());

	std::cout << "   " << solver_control.last_step()
			<< " CG iterations needed to obtain convergence."
			<< std::endl;
}

// Private method: output_results
template<int dim>
void ThermalProblem<dim>::output_results() const
{
	DataOut<dim> data_out;

	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, "solution");

	data_out.build_patches ();

	std::ofstream output (dim == 2 ?
			"solution-2d.vtk" :
			"solution-3d.vtk");

	data_out.write_vtk (output);
}


#endif /* THERMALPROBLEM_H_ */
