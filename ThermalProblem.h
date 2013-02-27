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

#include "BoundaryValues.h"

using namespace dealii;

template <int dim>
class ThermalProblem
{

public:
	ThermalProblem () {};
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


#endif /* THERMALPROBLEM_H_ */
