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
#include <deal.II/lac/block_sparse_matrix.h>
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
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

// General classes
#include "BoundaryGeometry.h"
#include "ScriptReader.h"
#include "TimeContainer.h"
#include "TMPostProcessor.h"
#include "Utility.h"
#include "VectorBoundary.h"

// Classes for Thermal
#include "ThermalRightHandSide.h"

// Classes for Elasticity
#include "Material.h"
#include "IsotropicMaterial.h"
#include "ElasticityRightHandSide.h"

#include <fstream>
#include <iostream>
#include <iterator>

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
	unsigned int n_comp;

	// Transient variables
	bool is_transient;
	TimeContainer tc;
	double initial_temperature;

	Triangulation<dim>  *triangulation;
	FESystem<dim>       fe;
	DoFHandler<dim>		dof_handler;
	ConstraintMatrix    hanging_node_constraints;

	std::vector<BoundaryGeometry<dim> *> * boundaries;
	std::vector<VectorBoundary<dim> *> temperature_bcs;
	std::vector<VectorBoundary<dim> *> flux_bcs;
	std::vector<VectorBoundary<dim> *> displacement_bcs;
	std::vector<VectorBoundary<dim> *> traction_bcs;

	std::vector< Material<dim> > materials;

	BlockSparsityPattern sparsity_pattern;
	BlockSparseMatrix<double> system_matrix;
	BlockVector<double> solution;
	BlockVector<double> prev_solution;
	BlockVector<double> system_rhs;
};

// Constructor
template<int dim>
ThermalElasticityProblem<dim>::ThermalElasticityProblem(Triangulation<dim> *triag) : fe(FE_Q<dim>(1), 1, FE_Q<dim>(1), dim), dof_handler(*triag) {
	triangulation = triag;
	verbosity = MIN_V;
	boundaries = 0;
	n_comp = fe.n_components();

	is_transient = false;
	initial_temperature = 0.0;
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
		Assert(tokens.size() == 6, ExcMessage("IsotropicMaterial command in the input script did not meet the expected parameters of material_id lambda mu alpha k."))
		// Expected tokens are:
		// material_id lambda mu alpha k
		materials.push_back(IsotropicMaterial<dim>(atoi(tokens[1].c_str()), atof(tokens[2].c_str()), atof(tokens[3].c_str()), atof(tokens[4].c_str()), atof(tokens[5].c_str())));
	}
	else if (tokens[0] == "SetVerbosity") {
		// There has already been an assertion in the Executive class
		verbosity = (Verbosity)atoi(tokens[1].c_str());
	}
	else if (tokens[0] == "TransientAnalysis") {
		Status("Analysis has been set to transient.", verbosity, MIN_V);
		Assert(tokens.size() == 4, ExcMessage("TransientAnalysis command in the input script did not meet the expected parameters of start_time end_time n_time_steps."))
		// Expected tokens are:
		// start_time end_time n_time_steps
		double start_time = atof(tokens[1].c_str());
		double end_time = atof(tokens[2].c_str());
		int n_time_steps = atoi(tokens[3].c_str());
		tc.initialize(start_time, end_time, n_time_steps);
		is_transient = true;
	}
	else if (tokens[0] =="InitialTemperature") {
		Assert(tokens.size() == 2, ExcMessage("InitialTemperature expects one arguement for the initial temperature."));
		// A single argument is for a constant temperature
		initial_temperature = atof(tokens[1].c_str());
	}
	else {
		return false;
	}

	return true;
} //process_cmd

// Public method: process_bc
template<int dim>
bool ThermalElasticityProblem<dim>::process_bc(std::vector<std::string> tokens)
{
	// Note this is a bit diconcerning, but the components and size of the finite element space must be know here so it will be hard coded
	// This must be changed if the FESystem is changed in its design
	// Right now thermal_comps = {0}
	//           elastic_comps = {1, ..., dim}
	std::vector<unsigned int> thermal_comps;
	thermal_comps.push_back(0);
	std::vector<unsigned int> elastic_comps;
	for (unsigned int i = 1; i < n_comp; i++)
		elastic_comps.push_back(i);


	if (tokens[0] == "DisplacementBoundary") {
		std::vector<double> disp_values;

		Assert(tokens.size() == 2 + dim, ExcMessage("The expected token number for DisplacementBoundary command was not met."))

		// Expected tokens for 2d are:
		// boundary_id x_value y_value
		// Expected tokens for 3d are:
		// boundary_id x_value y_value z_value
		for (int i = 0; i < dim; i++)
			disp_values.push_back(atof(tokens[2 + i].c_str()));

		VectorBoundary<dim> * bc_p = new VectorBoundary<dim>(atoi(tokens[1].c_str()), disp_values, elastic_comps, n_comp);
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

		VectorBoundary<dim> * bc_p = new VectorBoundary<dim>(atoi(tokens[1].c_str()),tract_values, elastic_comps, n_comp);
		traction_bcs.push_back(bc_p);
	}
	else if (tokens[0] == "TemperatureBoundary") {
		// Expected tokens for 2d and 3d are:
		// boundary_id value
		std::vector<double> temp_values;
		temp_values.push_back(atof(tokens[2].c_str()));

		VectorBoundary<dim> * bc_p = new VectorBoundary<dim>(atoi(tokens[1].c_str()), temp_values, thermal_comps, n_comp);
		temperature_bcs.push_back(bc_p);
	}
	else if (tokens[0] == "FluxBoundary") {
		// Expected tokens for 2d and 3d are:
		// boundary_id value
		std::vector<double> flux_values;
		flux_values.push_back(atof(tokens[2].c_str()));

		VectorBoundary<dim> * bc_p = new VectorBoundary<dim>(atoi(tokens[1].c_str()), flux_values, thermal_comps, n_comp);
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

	if (is_transient)
		std::cout << "\n\nAnalyzing time: " << tc.get_current_time() << "\n";

	assemble_system ();
	solve ();
	output_results ();

	if (is_transient) {
		while (tc.increment_time()) {
			std::cout << "\n\nAnalyzing time: " << tc.get_current_time() << "\n";
			assemble_system ();
			solve ();
			output_results ();
		}

		std::cout << "\nFinished transient analysis.\n";
	}
} // run

// Private method: setup_system
template<int dim>
void ThermalElasticityProblem<dim>::setup_system() // TODO: document this section
{
	Status("Starting setup_system.", verbosity, MIN_V);

	dof_handler.distribute_dofs (fe);
	DoFRenumbering::component_wise (dof_handler);

	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << "\n";

	std::vector<unsigned int> dofs_per_component (dim+1);
	DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);

	// Output for debug purposes
	if (verbosity >= DEBUG_V) {
		std::cout << "Dof's Per Component:" << "\n";
		for (unsigned int i = 0; i < dofs_per_component.size(); i++) {
			std::cout << "[" << i << "] " << dofs_per_component[i] << "\n";
		}
	}

	// Loop over indicies 1..dim+1 to get the dof for the displacements
	unsigned int n_u_nc = 0;
	for (unsigned int i = 1; i < dim + 1; i++)
		n_u_nc += dofs_per_component[i];
	// Number of dof for temperature (t) and displacement (u)
	const unsigned int n_t = dofs_per_component[0], n_u = n_u_nc;
	const unsigned int n_couplings = dof_handler.max_couplings_between_dofs();

	std::cout << n_t << " dims for temperature and " << n_u << " dims for displacement" << '\n';

	// Deal with the hanging nodes
	//hanging_node_constraints.clear();
	//DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
	//hanging_node_constraints.close();

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
	//hanging_node_constraints.condense (sparsity_pattern);
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

	if (is_transient) {
		// Create previous temperature vector
		prev_solution.reinit (2);
		prev_solution.block(0).reinit (n_t);
		prev_solution.block(1).reinit (n_u);
		prev_solution.collect_sizes ();

		// Assign initial temperature
		if (tc.index() == 0) {
			for (unsigned int i = 0; i < solution.block(0).size(); i++)
				solution.block(0)[i] = initial_temperature;
		}

	}

	Status("Completed setup_system.", verbosity, MIN_V);
} // setup_system

// Private method: assemble_system
template<int dim>
void ThermalElasticityProblem<dim>::assemble_system()
{
	Status("Starting assemble_system.", verbosity, MIN_V);

	Status("Starting initialization of assembly variables.", verbosity, MAX_V);
	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);

	const FEValuesExtractors::Scalar t_extract(0);
	const FEValuesExtractors::Vector u_extract(1);

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
	std::vector<double> local_prev_solution (dofs_per_cell);

	const ThermalRightHandSide<dim> thermal_rhs;
	const ElasticityRightHandSide<dim> elastic_rhs;
	std::vector<Vector<double> > elastic_rhs_values (n_q_points, Vector<double>(dim));
	Status("Completed initialization of assembly variables.", verbosity, MAX_V);

	Status("Starting the cell loop in assembly.", verbosity, MAX_V);
	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
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
				}

		fe_values.reinit (cell);
		cell_matrix = 0;
		cell_rhs = 0;
		cell->get_dof_indices (local_dof_indices);

		// Copy the right hand side terms for the elasticity to a vector to be used in this cell's assembly
		elastic_rhs.vector_value_list (fe_values.get_quadrature_points(), elastic_rhs_values);

		// If transient get the previous solution
		if (is_transient) {
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				local_prev_solution[i] = prev_solution(local_dof_indices[i]);
		}

		// Get cell's material properties
		SymmetricTensor<4, dim> stiffness;
		SymmetricTensor<2, dim> alpha;
		Tensor<2, dim> k;

		bool mat_found = false;
		for (unsigned int mat_ind = 0; mat_ind < materials.size(); mat_ind++)
			if (materials[mat_ind].get_id() == cell->material_id()) {
				stiffness = materials[mat_ind].get_stiffness();
				alpha = materials[mat_ind].get_expansion();
				k = materials[mat_ind].get_conduction();
				mat_found = true;
				break;
			}
		Assert(mat_found, ExcMessage("Material not found in assembly."))

		// Loop over the quad points of the cell and assemble the matrix and rhs
		// See http://www.dealii.org/developer/doxygen/deal.II/group__vector__valued.html end of section "An alternative approach"
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		{
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				// Temperature shape function for i
				const double phi_i_t = fe_values[t_extract].value(i, q_point);
				const Tensor<1,dim>  grad_phi_i_t = fe_values[t_extract].gradient(i, q_point);
				// Displacement shape function for i
				const Tensor<1,dim>  phi_i_u = fe_values[u_extract].value(i, q_point);
				const SymmetricTensor<2,dim>  grad_phi_i_u = fe_values[u_extract].symmetric_gradient(i, q_point);
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					// Temperature shape function for j
					const double phi_j_t = fe_values[t_extract].value(j, q_point);
					const Tensor<1,dim>  grad_phi_j_t = fe_values[t_extract].gradient(j, q_point);
					// Displacement shape function for j
					const SymmetricTensor<2,dim>  grad_phi_j_u = fe_values[u_extract].symmetric_gradient(j, q_point);

					cell_matrix(i,j) += (grad_phi_i_t * k * grad_phi_j_t
							+ -1*(grad_phi_i_u * (stiffness * alpha) * phi_j_t)
							+ grad_phi_i_u * (stiffness * grad_phi_j_u)
							) *
							fe_values.JxW (q_point);

					// Add transient term if this is not the first time step, since there is no previous time step
					if (is_transient  && tc.index() > 0) {
						cell_matrix(i,j) += (1/(tc.get_current_time() - tc.get_prev_time())) * phi_i_t * phi_j_t * fe_values.JxW (q_point);
						cell_rhs(i) += (1/(tc.get_current_time() - tc.get_prev_time())) * phi_i_t * phi_j_t * local_prev_solution[j] * fe_values.JxW (q_point);
				}

				cell_rhs(i) += (phi_i_t * thermal_rhs.value (fe_values.quadrature_point (q_point))
							+ phi_i_u * elastic_rhs_values[q_point]
							) * fe_values.JxW (q_point);
				}
			}
		}

		// Loop over faces and apply Nuemman BC's
		for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {
			if (cell->face(face)->at_boundary()) {
				// Apply Nuemman BC for thermal flux
				for (unsigned int flux_bc_index = 0; flux_bc_index < flux_bcs.size(); flux_bc_index++) {
					if (cell->face(face)->boundary_indicator() == flux_bcs[flux_bc_index]->get_id()) {
						fe_face_values.reinit(cell, face);
						for (unsigned int q_point = 0; q_point<n_face_q_points; ++q_point) {
							double flux_value = flux_bcs[flux_bc_index]->value(fe_face_values.quadrature_point(q_point));

							for (unsigned int i=0; i<dofs_per_cell; ++i) {
								const double phi_i_t = fe_values[t_extract].value(i, q_point);

								cell_rhs(i) += phi_i_t * flux_value * fe_face_values.JxW(q_point);
							}
						}
					}
				}

				// Apply Nuemman BC for traction
				for (unsigned int traction_bc_index = 0; traction_bc_index < traction_bcs.size(); traction_bc_index++) {
					if (cell->face(face)->boundary_indicator() == traction_bcs[traction_bc_index]->get_id()) {
						fe_face_values.reinit(cell, face);
						for (unsigned int q_point = 0; q_point<n_face_q_points; ++q_point) {
							Vector<double> traction_value(dim);
							// Be aware that if the fe formulation changes, this will need to be updated to the correct indicies
							traction_bcs[traction_bc_index]->vector_value(fe_face_values.quadrature_point(q_point), traction_value, 1, dim);

							for (unsigned int i=0; i<dofs_per_cell; ++i) {
								const Tensor<1,dim>  phi_i_u = fe_values[u_extract].value(i, q_point);
								double tmp_sum = 0;
								for (unsigned int j=0; j < dim; j++)
									tmp_sum += phi_i_u[j] * traction_value[j];
								cell_rhs(i) += tmp_sum * fe_face_values.JxW(q_point);
							}
						}
					}
				}
			}
		} // looping over faces

		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			for (unsigned int j=0; j<dofs_per_cell; ++j)
				system_matrix.add (local_dof_indices[i],
							local_dof_indices[j],
							cell_matrix(i,j));

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}
	//Status("Completed the cell loop in assembly.", verbosity, MAX_V);

	//hanging_node_constraints.condense(system_matrix);
	//hanging_node_constraints.condense(system_rhs);

	// TODO: Consider deriving rhs from TensorFunction rather than just Function, since it is more logical to describe as tensors
	// Apply Dirchlet BC for temperature
	ComponentMask temperature_mask = fe.component_mask(t_extract);
	for (unsigned int i = 0; i < temperature_bcs.size(); i++) {
		std::map<unsigned int,double> boundary_values_map;
		VectorTools::interpolate_boundary_values (dof_handler,
					temperature_bcs[i]->get_id(),
					*temperature_bcs[i],
					boundary_values_map,
					temperature_mask);

		MatrixTools::apply_boundary_values (boundary_values_map,
			  system_matrix,
			  solution,
			  system_rhs);
	}

	// Apply Dirchlet BC for displacement
	ComponentMask displacement_mask = fe.component_mask(u_extract);
	for (unsigned int i = 0; i < displacement_bcs.size(); i++) {
		std::map<unsigned int,double> boundary_values_map;
		VectorTools::interpolate_boundary_values (dof_handler,
					displacement_bcs[i]->get_id(),
					*displacement_bcs[i],
					boundary_values_map,
					displacement_mask);

		std::cout << "displacement_bc " << i << " " << *displacement_bcs[i] << "\n";

		MatrixTools::apply_boundary_values (boundary_values_map,
			  system_matrix,
			  solution,
			  system_rhs);
	}

	Status("Completed assemble_system.", verbosity, MIN_V);
} // assemble_system

template <int dim>
void ThermalElasticityProblem<dim>::solve ()
{
	Status("Starting solve for temperature.", verbosity, MIN_V);

	SolverControl           solver_control (1000, 1e-14);
	SolverCG<>              cg (solver_control);

	// If the problem is transient, the initial temperature is set and will be taken as the solution at time step 0
	if (!is_transient || tc.index() > 0) {
		// First solve for the temperatures
		cg.solve(system_matrix.block(0, 0), solution.block(0), system_rhs.block(0), PreconditionIdentity());
		//hanging_node_constraints.distribute(solution.block(0));
	}

	// Now solve for the displacements
	Vector<double> tmp (solution.block(1).size());
	system_matrix.block(1, 0).vmult(tmp, solution.block(0));
	system_rhs.block(1) -= tmp;
	tmp.reinit(0);

	Status("Starting solve for displacement.", verbosity, MIN_V);

	PreconditionSSOR<> preconditioner;
	preconditioner.initialize(system_matrix.block(1, 1), 1.2);
	cg.solve(system_matrix.block(1, 1), solution.block(1), system_rhs.block(1), preconditioner);
	//hanging_node_constraints.distribute(solution.block(1));

	std::cout << "Displacement:\n";
	for(unsigned int i = 0; i < solution.block(1).size(); i++)
		std::cout << i << ": " << solution.block(1)[i] << "\n";

	//Status("Completed solve for displacement.", verbosity, MIN_V);

	if (is_transient) {
		Status("Saving current solution for next time step.", verbosity, MIN_V);
		prev_solution = solution;
	}
} // solve

template <int dim>
void ThermalElasticityProblem<dim>::output_results () const
{
	TEPostProcessor<dim> post_process;
	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, post_process);
	data_out.build_patches();
	if (is_transient) {
		std::ofstream output ("solution" + tc.index_str() + ".vtk");
		data_out.write_vtk(output);
	}
	else {
		std::ofstream output ("solution.vtk");
		data_out.write_vtk(output);
	}
}


}

#endif /* THERMALELASTICITYPROBLEM_H_ */
