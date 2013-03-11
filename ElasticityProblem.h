/*
 * ElasticityProblem.h
 *
 *  Created on: Mar 7, 2013
 *      Author: kballard
 */

#ifndef ELASTICITYPROBLEM_H_
#define ELASTICITYPROBLEM_H_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include "BoundaryCondition.h"
#include "BoundaryGeometry.h"
#include "RightHandSide.h"
#include "ScriptReader.h"
#include "Utility.h"

#include <fstream>
#include <iostream>

namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class ElasticityProblem
{

public:
	ElasticityProblem ();
	~ElasticityProblem();

	void run_test ();
	void run (ScriptReader *script_reader);

private:
	void make_grid_test ();

	void setup_system();
	void assemble_system ();
	void solve ();
	void output_results () const;

	ScriptReader* sr;
	Verbosity verbosity;

	Triangulation<dim>  triangulation;
	FE_Q<dim>          	fe;
	DoFHandler<dim>		dof_handler;

	std::vector<BoundaryGeometry<dim> *> boundaries;
	std::vector<TemperatureBoundary<dim> *> temperature_bcs;
	std::vector<FluxBoundary<dim> *> flux_bcs;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double>       solution;
	Vector<double>       system_rhs;

	// State variables
	bool mesh_initialized;
};

}

#endif /* ELASTICITYPROBLEM_H_ */
