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
#include <deal.II/base/exceptions.h>
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
#include <string>
#include <vector>

#include "BoundaryCondition.h"
#include "BoundaryGeometry.h"
#include "ThermalRightHandSide.h"
#include "ScriptReader.h"

#include "TemperatureBoundary.h"
#include "FluxBoundary.h"

#include "Utility.h"


namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class ThermalProblem
{

public:
	ThermalProblem (Triangulation<dim> *triag);
	~ThermalProblem();

	bool process_cmd(std::vector<std::string> tokens);
	bool process_bc(std::vector<std::string> tokens);

	void run(std::vector<BoundaryGeometry<dim> *> *bound);

	Vector<double>* get_thermal_sol();
	std::vector<std::vector<Tensor< 1, dim>>>* get_thermal_grad();

private:
	void setup_system();
	void assemble_system ();
	void solve ();
	void output_results () const;

	Verbosity verbosity;

	Triangulation<dim>  *triangulation;
	FE_Q<dim>          	fe;
	DoFHandler<dim>		dof_handler;

	std::vector<BoundaryGeometry<dim> *> * boundaries;
	std::vector<TemperatureBoundary<dim> *> temperature_bcs;
	std::vector<FluxBoundary<dim> *> flux_bcs;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double>       solution;
	Vector<double>       system_rhs;
	std::vector<std::vector<Tensor< 1, dim>>> thermal_grads;

};

// Constructor
template<int dim>
ThermalProblem<dim>::ThermalProblem(Triangulation<dim> *triag) : fe(1), dof_handler(*triag) {
	triangulation = triag;
	verbosity = MIN_V;
	boundaries = 0;
}

// Destructor
template<int dim>
ThermalProblem<dim>::~ThermalProblem() {
	for (unsigned int i = 0; i < temperature_bcs.size(); i++) {
		delete temperature_bcs[i];
	}

	for (unsigned int i = 0; i < flux_bcs.size(); i++) {
		delete flux_bcs[i];
	}

	dof_handler.clear ();
}

// Public method: process_cmd
template<int dim>
bool ThermalProblem<dim>::process_cmd(std::vector<std::string> tokens)
{
	// No commands right now
	return false;
}

// Public method: process_bc
template<int dim>
bool ThermalProblem<dim>::process_bc(std::vector<std::string> tokens)
{
	if (tokens[0] == "TemperatureBoundary") {
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
}

// Public method: get_thermal_sol
template<int dim>
Vector<double>* ThermalProblem<dim>::get_thermal_sol()
{
	return &solution;
}

// Public method: get_thermal_grad
template<int dim>
std::vector<std::vector<Tensor< 1, dim>>>* ThermalProblem<dim>::get_thermal_grad()
{
	return &thermal_grads;
}

// Public method: run
template<int dim>
void ThermalProblem<dim>::run(std::vector<BoundaryGeometry<dim> *> *bound)
{
	std::cout << "Solving test problem in " << dim << " space dimensions." << std::endl;

	boundaries = bound;

	setup_system ();
	assemble_system ();
	output_results ();
}

// Private method: setup_system
template<int dim>
void ThermalProblem<dim>::setup_system()
{
	Status("Starting setup_system.", verbosity, MIN_V);

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
	Status("Starting assemble_system.", verbosity, MIN_V);

	Status("Starting initialization of assembly variables.", verbosity, MAX_V);
	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);

	const ThermalRightHandSide<dim> right_hand_side;

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

		// Loop over faces and apply Nuemman BC's
		for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {
			if (cell->face(face)->at_boundary()) {
				for (unsigned int flux_bc_index = 0; flux_bc_index < flux_bcs.size(); flux_bc_index++) {
					if (cell->face(face)->boundary_indicator() == flux_bcs[flux_bc_index]->get_id()) {
						fe_face_values.reinit(cell, face);
						for (unsigned int q_point = 0; q_point<n_face_q_points; ++q_point) {
							const double neumann_value = flux_bcs[flux_bc_index]->value(fe_face_values.quadrature_point(q_point));

							for (unsigned int i=0; i<dofs_per_cell; ++i)
								cell_rhs(i) += (neumann_value * fe_face_values.shape_value(i,q_point) * fe_face_values.JxW(q_point));
						}
					}
				}
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
	Status("Completed the cell loop in assembly.", verbosity, MAX_V);


	for (unsigned int i = 0; i < temperature_bcs.size(); i++) {
		std::map<unsigned int,double> boundary_values_map;
		VectorTools::interpolate_boundary_values (dof_handler,
					temperature_bcs[i]->get_id(),
					*temperature_bcs[i],
					boundary_values_map);

		MatrixTools::apply_boundary_values (boundary_values_map,
			  system_matrix,
			  solution,
			  system_rhs);
	}

	// Call the solve
	solve();

	// Compute the thermal_solution gradients
	cell = dof_handler.begin_active();
	for (; cell!=endc; ++cell)
	{
		fe_values.reinit (cell);

		std::vector<Tensor< 1, dim>> cell_grads;
		cell_grads.resize(n_q_points);

		fe_values.get_function_gradients(solution, cell_grads);

		thermal_grads.push_back(cell_grads);
	}
}

// Private method: solve
template<int dim>
void ThermalProblem<dim>::solve()
{
	Status("Starting solve.", verbosity, MIN_V);

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
	Status("Starting output_results.", verbosity, MIN_V);

	DataOut<dim> data_out;

	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, "solution");

	data_out.build_patches ();

	std::ofstream output (dim == 2 ?
			"thermal-solution-2d.vtk" :
			"thermal-solution-3d.vtk");

	data_out.write_vtk (output);
}

}

#endif /* THERMALPROBLEM_H_ */
