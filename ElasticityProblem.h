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
class ElasticityProblem
{

public:
	ElasticityProblem (Triangulation<dim> *triag);
	~ElasticityProblem();

	bool process_cmd(std::vector<std::string> tokens);
	bool process_bc(std::vector<std::string> tokens);

	Vector<double>& run(std::vector<BoundaryGeometry<dim> *> *bound);

private:
	void setup_system();
	void assemble_system ();
	void solve ();
	void output_results () const;

	Verbosity verbosity;

	Triangulation<dim>  *triangulation;
	FE_Q<dim>          	fe;
	DoFHandler<dim>		dof_handler;
	ConstraintMatrix     hanging_node_constraints;

	std::vector<BoundaryGeometry<dim> *> * boundaries;
	std::vector<DisplacementBoundary<dim> *> displacement_bcs;
	std::vector<TractionBoundary<dim> *> traction_bcs;

	std::vector<IsotropicMaterial> materials;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double>       solution;
	Vector<double>       system_rhs;
};

// Constructor
template<int dim>
ElasticityProblem<dim>::ElasticityProblem(Triangulation<dim> *triag) : fe(1), dof_handler(*triag) {
	triangulation = triag;
	verbosity = MAX_V;
	boundaries = 0;
}

// Destructor
template<int dim>
ElasticityProblem<dim>::~ElasticityProblem() {
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
bool ElasticityProblem<dim>::process_cmd(std::vector<std::string> tokens)
{
	// TODO: Add support for orthotropic materials
	if (tokens[0] == "IsotropicMaterial") {
		Status("Adding an isotropic material.", verbosity, MAX_V);
		// Expected tokens are:
		// material_id lambda mu
		materials.push_back(IsotropicMaterial(atoi(tokens[1].c_str()), atof(tokens[2].c_str()), atof(tokens[3].c_str())));
	}
	else {
		return true;
	}

	return true;
} //process_cmd

// Public method: process_bc
template<int dim>
bool ElasticityProblem<dim>::process_bc(std::vector<std::string> tokens)
{
	if (tokens[0] == "DisplacementBoundary") {
		std::vector<double> disp_values;

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

		// Expected tokens for 2d are:
		// boundary_id x_value y_value
		// Expected tokens for 3d are:
		// boundary_id x_value y_value z_value
		for (int i = 0; i < dim; i++)
			tract_values.push_back(atof(tokens[2 + i].c_str()));

		TractionBoundary<dim> * bc_p = new TractionBoundary<dim>(atoi(tokens[1].c_str()),tract_values);
		traction_bcs.push_back(bc_p);
	}
	else {
		return false;
	}

	return true;
} // process_bc

// Public method: run
template<int dim>
Vector<double>& ElasticityProblem<dim>::run(std::vector<BoundaryGeometry<dim> *> *bound)
{
	std::cout << "Solving test problem in " << dim << " space dimensions." << std::endl;

	boundaries = bound;

	setup_system ();
	assemble_system ();
	solve ();
	output_results ();

	return solution;
} // run

// Private method: setup_system
template<int dim>
void ElasticityProblem<dim>::setup_system() // TODO: document this section
{
	Status("Starting setup_system.", verbosity, MIN_V);

	dof_handler.distribute_dofs (fe);

	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << "\n";

	// Deal with the hanging nodes
	hanging_node_constraints.clear();
	DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
	hanging_node_constraints.close();

	sparsity_pattern.reinit (dof_handler.n_dofs(),
							 dof_handler.n_dofs(),
							 dof_handler.max_couplings_between_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

	hanging_node_constraints.condense (sparsity_pattern);

	sparsity_pattern.compress();

	system_matrix.reinit (sparsity_pattern);

	solution.reinit (dof_handler.n_dofs());
	system_rhs.reinit (dof_handler.n_dofs());
} // setup_system

// Private method: assemble_system
template<int dim>
void ElasticityProblem<dim>::assemble_system()
{
	Status("Starting assemble_system.", verbosity, MIN_V);

	Status("Starting initialization of assembly variables.", verbosity, MAX_V);
	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);

	const ElasticityRightHandSide<dim> right_hand_side;

	FEValues<dim> fe_values (fe, quadrature_formula,
			update_values   | update_gradients |
			update_quadrature_points | update_JxW_values);

	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
			update_values         | update_quadrature_points  |
			update_normal_vectors | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	// Used to modify the right hand side for the constraints
	std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim));
	Status("Completed initialization of assembly variables.", verbosity, MAX_V);

	Status("Starting the cell loop in assembly.", verbosity, MAX_V);
	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
				endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		// Apply boundary_indicators to the faces given the boundary conditions
		// Check if the center of the face is on one of the given boundaries and set its boundary indicator
		// In the case it is not, the boundary indicator will remain 0
		for (unsigned int face = 0; face<GeometryInfo<dim>::faces_per_cell; face++)
			for (unsigned int boundary_index = 0; boundary_index < boundaries->size(); boundary_index++)
				if ((*boundaries)[boundary_index]->point_on_geometry(cell->face(face)->center())) {
					cell->face(face)->set_boundary_indicator((*boundaries)[boundary_index]->get_id()); // TODO: Should this be set_all_boundary_indicators?
					if (verbosity == DEBUG_V) {
						std::cout << "Setting boundary id " << (*boundaries)[boundary_index]->get_id() << " at face " << face << " with center: (";
						for (int i = 0; i < dim; i++)
							std::cout << cell->face(face)->center()(i) << ", ";
						std::cout << ")" << std::endl;
					}
				}

		fe_values.reinit (cell);
		cell_matrix = 0;
		cell_rhs = 0;

		right_hand_side.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
		// Assemble stiffness matrix
		for (unsigned int i=0; i<dofs_per_cell; i++)
		{
			const unsigned int component_i = fe.system_to_component_index(i).first;
			for (unsigned int j=0; j<dofs_per_cell; j++)
			{
				const unsigned int component_j = fe.system_to_component_index(j).first;
				for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
					double lambda = 0;
					double mu = 0;

					bool mat_found = false;
					for (unsigned int mat_ind = 0; mat_ind < materials.size(); mat_ind++)
						if (materials[mat_ind].get_id() == cell->material_id()) {
							lambda = materials[mat_ind].get_lambda();
							mu = materials[mat_ind].get_mu();
							mat_found = true;
							break;
						}

					Assert(mat_found, ExcMessage("Material not found in assembly."))


					cell_matrix(i,j) +=(
	                        (fe_values.shape_grad(i,q_point)[component_i] * fe_values.shape_grad(j,q_point)[component_j] * lambda)
	                        + (fe_values.shape_grad(i,q_point)[component_j] * fe_values.shape_grad(j,q_point)[component_i] * mu)
	                        + ((component_i == component_j) ? (fe_values.shape_grad(i,q_point) * fe_values.shape_grad(j,q_point) * mu)  : 0)
	                      ) * fe_values.JxW(q_point);
				}
			}
		}

		// Assemble the rhs
		for (unsigned int i=0; i<dofs_per_cell; i++)
		{
			const unsigned int component_i = fe.system_to_component_index(i).first;
			for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
				cell_rhs(i) += (fe_values.shape_value (i, q_point) *
						rhs_values[q_point](component_i) *
						fe_values.JxW (q_point));
			}
		}

		// Loop over faces and apply Nuemman BC's
		for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {
			if (cell->face(face)->at_boundary()) {
				for (unsigned int traction_bc_index = 0; traction_bc_index < traction_bcs.size(); traction_bc_index++) {
					if (cell->face(face)->boundary_indicator() == traction_bcs[traction_bc_index]->get_id()) {
						fe_face_values.reinit(cell, face);
						for (unsigned int q_point = 0; q_point<n_face_q_points; ++q_point) {
							const double neumann_value = traction_bcs[traction_bc_index]->value(fe_face_values.quadrature_point(q_point));

							for (unsigned int i=0; i<dofs_per_cell; ++i)
								cell_rhs(i) += (neumann_value * fe_face_values.shape_value(i,q_point) * fe_face_values.JxW(q_point));
						}
					}
				}
			}
		}

		// Map values to global matrix and rhs
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
	Status("Completed the cell loop in assembly.", verbosity, MAX_V);


	for (unsigned int i = 0; i < displacement_bcs.size(); i++) {
		std::map<unsigned int,double> boundary_values_map;
		VectorTools::interpolate_boundary_values (dof_handler,
					displacement_bcs[i]->get_id(),
					*displacement_bcs[i],
					boundary_values_map);

		MatrixTools::apply_boundary_values (boundary_values_map,
			  system_matrix,
			  solution,
			  system_rhs);
	}
} // assemble_system

template <int dim>
void ElasticityProblem<dim>::solve ()
{
	SolverControl           solver_control (1000, 1e-12);
	SolverCG<>              cg (solver_control);

	PreconditionSSOR<> preconditioner;
	preconditioner.initialize(system_matrix, 1.2);

	cg.solve (system_matrix, solution, system_rhs,
			  preconditioner);

	hanging_node_constraints.distribute (solution);
}

template <int dim>
void ElasticityProblem<dim>::output_results () const
{
	std::string filename = "elastic-solution.vtk";
	std::ofstream output (filename.c_str());

	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);



	std::vector<std::string> solution_names;
	switch (dim)
	{
	case 1:
		solution_names.push_back ("displacement");
		break;
	case 2:
		solution_names.push_back ("x_displacement");
		solution_names.push_back ("y_displacement");
		break;
	case 3:
		solution_names.push_back ("x_displacement");
		solution_names.push_back ("y_displacement");
		solution_names.push_back ("z_displacement");
		break;
	default:
		Assert (false, ExcNotImplemented());
		break;
	}

	data_out.add_data_vector (solution, solution_names);
	data_out.build_patches ();
	data_out.write_vtk (output);
}

} // namespace

#endif /* ELASTICITYPROBLEM_H_ */
