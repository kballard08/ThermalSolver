/*
 * ThermalElasticityProblem.h
 *
 *  Created on: Mar 26, 2013
 *      Author: kballard
 */

#ifndef THERMALELASTICITYPROBLEM_H_
#define THERMALELASTICITYPROBLEM_H_

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
#include "DisplacementBoundary.h"
#include "TractionBoundary.h"
#include "ElasticityRightHandSide.h"
#include "IsotropicMaterial.h"
#include "ScriptReader.h"
#include "Utility.h"

#include <fstream>
#include <iostream>

namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class ThermalElasticityProblem
{
public:
	ThermalElasticityProblem(Triangulation<dim> *triag);
	~ThermalElasticityProblem();

	bool process_cmd(std::vector<std::string> tokens);
	bool process_bc(std::vector<std::string> tokens);

	void run(std::vector<BoundaryGeometry<dim> *> *bound);

private:
	void setup_system();
	void assemble_system ();
	void solve ();
	void output_results () const;

	Verbosity verbosity;

	Triangulation<dim>  *triangulation;
	FESystem<dim>       fe;
	DoFHandler<dim>		dof_handler;
	ConstraintMatrix    hanging_node_constraints;

	std::vector<BoundaryGeometry<dim> *> * boundaries;
	std::vector<TemperatureBoundary<dim> *> temperature_bcs;
	std::vector<FluxBoundary<dim> *> flux_bcs;
	std::vector<DisplacementBoundary<dim> *> displacement_bcs;
	std::vector<TractionBoundary<dim> *> traction_bcs;

	std::vector<IsotropicMaterial> materials;

	BlockSparsityPattern sparsity_pattern;
	BlockSparseMatrix<double> system_matrix;
	BlockVector<double> solution;
	BlockVector<double> system_rhs;
};

// Constructor
template<int dim>
ThermalElasticityProblem<dim>::ThermalElasticityProblem(Triangulation<dim> *triag) : fe(FE_Q<dim>(1), 1, FE_Q<dim>(1), dim), dof_handler(*triag) {
	triangulation = triag;
	verbosity = MAX_V;
	boundaries = 0;
}

// Destructor
template<int dim>
ThermalElasticityProblem<dim>::~ThermalElasticityProblem() {
	for (unsigned int i = 0; i < temperature_bcs.size(); i++) {
		delete temperature_bcs[i];
	}

	for (unsigned int i = 0; i < flux_bcs.size(); i++) {
		delete flux_bcs[i];
	}

	for (unsigned int i = 0; i < displacement_bcs.size(); i++) {
		delete displacement_bcs[i];
	}

	for (unsigned int i = 0; i < traction_bcs.size(); i++) {
		delete traction_bcs[i];
	}

	dof_handler.clear ();
}

// Public method: process_cmd
template<int dim>
bool ThermalElasticityProblem<dim>::process_cmd(std::vector<std::string> tokens)
{
	// TODO: Add support for orthotropic materials
	if (tokens[0] == "IsotropicMaterial") {
		Status("Adding an isotropic material.", verbosity, MAX_V);
		Assert(tokens.size() == 5, ExcMessage("IsotropicMaterial command in the input script did not meet the expected parameters of material_id lambda mu alpha."))
		// Expected tokens are:
		// material_id lambda mu alpha
		materials.push_back(IsotropicMaterial(atoi(tokens[1].c_str()), atof(tokens[2].c_str()), atof(tokens[3].c_str()), atof(tokens[4].c_str())));
	}
	else {
		return true;
	}

	return true;
} //process_cmd

// Public method: process_bc
template<int dim>
bool ThermalElasticityProblem<dim>::process_bc(std::vector<std::string> tokens)
{
	if (tokens[0] == "DisplacementBoundary") {
		std::vector<double> disp_values;

		Assert(tokens.size() == 2 + dim, ExcMessage("The expected token number for DisplacementBoundary command was not met."))

		// Expected tokens for 2d are:
		// boundary_id x_value y_value
		// Expected tokens for 3d are:
		// boundary_id x_value y_value z_value
		for (int i = 0; i < dim; i++)
			disp_values.push_back(atof(tokens[2 + i].c_str()));

		DisplacementBoundary<dim> * bc_p = new DisplacementBoundary<dim>(atoi(tokens[1].c_str()), disp_values);
		displacement_bcs.push_back(bc_p);
	}
	else if (tokens[0] == "TractionBoundary") {
		std::vector<double> tract_values;

		Assert(tokens.size() == 2 + dim, ExcMessage("The expected token number for TractionBoundary command was not met."))

		// Expected tokens for 2d are:
		// boundary_id x_value y_value
		// Expected tokens for 3d are:
		// boundary_id x_value y_value z_value
		for (int i = 0; i < dim; i++)
			tract_values.push_back(atof(tokens[2 + i].c_str()));

		TractionBoundary<dim> * bc_p = new TractionBoundary<dim>(atoi(tokens[1].c_str()),tract_values);
		traction_bcs.push_back(bc_p);
	}
	else if (tokens[0] == "TemperatureBoundary") {
		// Expected tokens for 2d and 3d are:
		// boundary_id value
		TemperatureBoundary<dim> * bc_p = new TemperatureBoundary<dim>(atoi(tokens[1].c_str()), atof(tokens[2].c_str()));
		temperature_bcs.push_back(bc_p);
	}
	else if (tokens[0] == "FluxBoundary") {
		// Expected tokens for 2d and 3d are:
		// boundary_id value
		FluxBoundary<dim> * bc_p = new FluxBoundary<dim>(atoi(tokens[1].c_str()), atof(tokens[2].c_str()));
		flux_bcs.push_back(bc_p);
	}
	else {
		return false;
	}

	return true;
} // process_bc

// Public method: run
template<int dim>
void ThermalElasticityProblem<dim>::run(std::vector<BoundaryGeometry<dim> *> *bound)
{
	std::cout << "Solving test problem in " << dim << " space dimensions." << std::endl;

	boundaries = bound;

	setup_system ();
	assemble_system ();
	solve ();
	output_results ();
} // run

// Private method: setup_system
template<int dim>
void ThermalElasticityProblem<dim>::setup_system() // TODO: document this section
{
	Status("Starting setup_system.", verbosity, MIN_V);

	dof_handler.distribute_dofs (fe);
	DoFRenumbering::component_wise (dof_handler);

	std::cout << "   Number of degrees o freedom: " << dof_handler.n_dofs() << "\n";

	std::vector<unsigned int> dofs_per_component (dim+1);
	DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
	// Number of dof for temperature (t) and displacement (u)
	// Is this neccessary?  I guess if different elements are used later
	const unsigned int n_t = dofs_per_component[0], n_u = dofs_per_component[1];
	const unsigned int n_couplings = dof_handler.max_couplings_between_dofs();

	// Deal with the hanging nodes
	hanging_node_constraints.clear();
	DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
	hanging_node_constraints.close();

	// Create sparsity pattern
	// [ thermal stiffness terms         , thermal terms from deformation ] * [ t ] = [ heat generation terms ]
	// [ displacement terms from heat    , displacement stiffness terms   ]   [ u ]   [ body force terms      ]
	sparsity_pattern.reinit (2,2);
	sparsity_pattern.block(0,0).reinit (n_t, n_t, n_couplings);
	sparsity_pattern.block(1,0).reinit (n_u, n_t, n_couplings);
	sparsity_pattern.block(0,1).reinit (n_t, n_u, n_couplings);
	sparsity_pattern.block(1,1).reinit (n_u, n_u, n_couplings);
	sparsity_pattern.collect_sizes();

	// Create and compress
	DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
	hanging_node_constraints.condense (sparsity_pattern);
	sparsity_pattern.compress();

	// Create solution and rhs vectors
	system_matrix.reinit (sparsity_pattern);
	solution.reinit (2);
	solution.block(0).reinit (n_t);
	solution.block(1).reinit (n_u);
	solution.collect_sizes ();
	system_rhs.reinit (2);
	system_rhs.block(0).reinit (n_t);
	system_rhs.block(1).reinit (n_u);
	system_rhs.collect_sizes ();

	Status("Completed setup_system.", verbosity, MIN_V);
} // setup_system

}

#endif /* THERMALELASTICITYPROBLEM_H_ */