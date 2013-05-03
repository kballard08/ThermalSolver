/*
 * ThermalElasticityProblem.h
 *
 *  Created on: Mar 26, 2013
 *      Author: kballard
 */

#ifndef THERMALELASTICITYPROBLEM_H_
#define THERMALELASTICITYPROBLEM_H_

/// Include files
///=============

/// Include applicable deal II headers
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

/// Include general classes
#include "BoundaryGeometry.h"
#include "ScriptReader.h"
#include "TimeContainer.h"
#include "TEDataOut.h"
#include "TEPostProcessor.h"
#include "Utility.h"
#include "VectorBoundary.h"

/// Classes neccessary for thermal problem
#include "ThermalRightHandSide.h"

/// Classes neccessary for elasticity problem
#include "Material.h"
#include "IsotropicMaterial.h"
#include "ElasticityRightHandSide.h"

/// Std headers
#include <fstream>
#include <iostream>
#include <iterator>

/// Declare the namespace for this program
namespace FEASolverNS
{

/// Use the deal ii namespace inside the declared namespace
using namespace dealii;

///ThermalElasticityProblem Class Description
///==========================================

/**
 * The ThermalElasticityProblem encapsulates the setup, assembly, solve, and output for the thermomechanical problem.
 * There are dependancies on other classes from included files for functionality, but the fe forumation and main
 * theory of the problem lies in this class.  The Executive instantiates this class and calls the run command.
 *
 * This class is based on step-4, step-8, and step-40
 */
template <int dim>
class ThermalElasticityProblem
{
public:
	/// Constructor and Destructor
	///--------------------------

	/**
	 * Constructor takes a pointer to the triangulation to be used for the problem.  The Executive class contains/owns
	 * the triangulation for reasons explained in the Executive class.  This class must not delete the pointer, as it does
	 * not own it.  The Executive class will handle the destruction of the triangulation object.
	 */
	ThermalElasticityProblem(Triangulation<dim> *triag);
	/**
	 * Destructor must clean up the memory used by this class.  Namely several vectors of pointers.  There are safer ways
	 * to handle this using the C++11 standard, but as long as care is taken storing a vector of pointers is fine since
	 * there is no chance for another class to delete it.
	 */
	~ThermalElasticityProblem();

	/// Public Methods
	///--------------

	/**
	 * There is sort of a remenant of the method I was used to for reading in files.  Deal II contains methods and structures
	 * for interfacing with a parameters file, but this would limit the program in one way.  Our research group has FEA analyzes
	 * that are driven by text files that have an established format.  The formats used in this program are not the same neccessarily
	 * but the concept and layout is very similiar.  This design choice will allow me to easily add support for reading in our
	 * existing FEA scripts and drive both code frameworks (Beta and this one).
	 */

	/**
	 * Reads commands that are relevant to the analysis, such as material information, initial temperature, verbosity of output,
	 * or options concerning the transient analysis.
	 *
	 * Input script syntax:
	 * IsotropicMaterial material_id E nu alpha k  (where the material_id corresponds to that of the mesh, E is Youngs's modulus,
	 * 		nu is Poisson's ratio, alpha is coefficient of thermal expansion, and k is coefficient of thermal conductivity)
	 * TransientAnalysis start end num_steps (where num_steps are the number of stepx excluding the first on as that is simply a
	 * 		calculation of only elasticity problem)
	 * InitialTemperature value  (Assumed to be 0 if not set by the script)
	 */
	bool process_cmd(std::vector<std::string> tokens);
	/**
	 * Reads boundary conditions for the analysis.  It assumes that every boundary id introduced will have a corresponding
	 * boundary condition.  (TODO: add assertion that this actually happens, since I cannot think of a reason that you would
	 * want otherwise)  Currently only boundary conditions for Dirchlet and Nuemman for both elasticity and thermal are supported
	 * Robin and other modified boundary condition types could be added here and implemented into the assembly methods.  The reason
	 * only the mentioned two were added was because in the analyzes my research group is conerned with, there is no neccessity for
	 * other types, but in the future they can be added as needed.
	 *
	 * Input script syntax:
	 * TemperatureBoundary boundary_id temperature_value
	 * FluxBoundary boundary_id flux_value
	 * DisplacementBoundary boundary_id x_value y_value (z_value for 3D)
	 * TractionBoundary boundary_id x_value y_value (z_value for 3D)
	 */
	bool process_bc(std::vector<std::string> tokens);

	/**
	 * The method called by the Executive class once the input script has been read and the problem is ready to be solved.
	 * The argument for a pointer boundaries comes from the Executive class.  See the Executive class for reasons why
	 * this is not inside the ThermalElastictyProblem class.
	 */
	void run(std::vector<BoundaryGeometry<dim> *> *bound);

private:
	/// Private Methods
	///---------------

	/**
	 * The setup_system method determines the number of dof's for the blocks of the system, creates sparsity patters,
	 * intializes the system vectors and matrices, initializes the temperature values, and applies boundary id's.
	 */
	void setup_system();
	/**
	 * The assemble_system forms the system matrix and right hand side vector, applies boundary conditions, and hanging node
	 * constraints. Note: in the case of a transient analysis, this method is only called onece by the run method.
	 * Subsequent assemble calls are done to the transient_assembly_system.
	 */
	void assemble_system ();
	/**
	 * The transient_assemble_system method only assembles the transient terms.  It sets the thermal block of the system matrix
	 * and right hand side vector to zero and recalculates the terms.  The elasticty blocks do not change across the time
	 * steps, so they are not recalculated to save time, which ended up being a significant time savings on the order of 20%
	 * for time steps after the intial assembly and solve.  The method also calls apply_dirchlet_bcs out of
	 * neccessity since the solution vector is modified during the solve.
	 */
	void transient_assemble_system ();
	/**
	 * The apply_dirchlet_bcs method does as its name implies.  It applies the dirchlet boundary conditions to the system matrix,
	 * right hand side vector, and solution vector.  This method is called by both assembly_system and transient_assemble_system.
	 */
	void apply_dirchlet_bcs ();
	/**
	 * The solve method solves the heat equation first, multiplies the thermal solution with the thermal-elasticity coupling block
	 * of the system matrix, subtracts the result from the right hand side and solves the elasticity problem.  The method also
	 * saves the previous solution for use at the next time step.
	 */
	void solve ();
	/**
	 * Outputs the results (temperature, displacement, strain, and stress) to a vtk file titled solution.  If the problem is
	 * transient, each file is followed by an index allowing the files to be read a single solution data source in Visit or
	 * Paraview.
	 */
	void output_results () const;

	/// Private Properties
	///------------------

	/// Simply a variable to keep track of how much to output to the screen.
	Verbosity verbosity;
	/**
	 * Total number of components in the fe system.  This is simply for ease during the private methods.  The variable is
	 * initialized in the constructor.
	 */
	unsigned int n_comp;

	/// Transient variables
	/// A boolean indicating if the analysis is transient for use in the assembly and output methods
	bool is_transient;
	/**
	 * I realized deal II had methods to keep track of time after the program had been designed.  The design was not changed
	 * since the TimeContainer class accomplished what is needed, and allows a few benefits.  See the TimeContainer class
	 * for more details.
	 */
	TimeContainer tc;
	/// The initial temperature as set by the input script (defaults to 0 if unset).  Is really only used in setup_system.
	double initial_temperature;

	/// A pointer to the Executive class's triangulation, see the constructor description and Executive class documentation
	/// for more information on this.
	Triangulation<dim>  *triangulation;
	/// The FESystem, which is block formulated
	FESystem<dim>       fe;
	/// The DofHandler, which handles the dof's attached to the triangulation.  This is owned by the class, unlike the
	/// triangulation.
	DoFHandler<dim>		dof_handler;
	/// ConstraintMatrix that contains the constraints due to hanging nodes.  In the future, the ConstraintMatrix can be
	/// extended to handle multipoint constraints, which are needed for periodic boundary conditions.
	ConstraintMatrix    hanging_node_constraints;

	/**
	 * Pointer to Executive class's vector of boundaries. This pointer is owned by the Executive class and should not
	 * be deleted here, as it will be destructed by the Executive class. Again, this should be implemented with a shared
	 * pointer later, but for now the code is careful.
	 */
	std::vector<BoundaryGeometry<dim> *> * boundaries;
	/**
	 * Vector of points to boundary conditions described by the VectorBoundary class.  These pointers are deleted in the
	 * destructor of this class.  Each boundary condition type are separted into their own vector.  These were designed
	 * to be vector of pointers when the program had separate classes for the elasticity and thermal problem.  Now the
	 * architecture could be changed to be a vector of VectorBoundary objects rather than pointers, which is better
	 * programming practice.
	 */
	/// Temperature boundary conditions
	std::vector<VectorBoundary<dim> *> temperature_bcs;
	/// Thermal flux boundary conditions
	std::vector<VectorBoundary<dim> *> flux_bcs;
	/// Displacement boundary conditions
	std::vector<VectorBoundary<dim> *> displacement_bcs;
	/// Traction (stress operating on the normal of the boundary) boundary conditions
	std::vector<VectorBoundary<dim> *> traction_bcs;

	/// A vector of Material objects. This is a general description of a material and can contain orthotropic to
	/// isotropic materials.
	std::vector< Material<dim> > materials;

	/// Sparsity pattern for the system
	BlockSparsityPattern sparsity_pattern;
	/// Block system matrix as described in the problems formulation.
	BlockSparseMatrix<double> system_matrix;
	/// The block solution vector containing the thermal and elasticity solutions.
	BlockVector<double> solution;
	/**
	 * The previous solution vector needed for transient terms. Only the previous temperature is needed, but during the
	 * assembly methods, the cell only knows the global indicies and doesn't know about which blocks they fall in.
	 * To allow the cell to get the previous solution corresponding to the dof's it is responsibile for, a full
	 * BlockVector is kept. This keeps it open for dof renumbering in the future without affecting the aseembly part
	 * of the code.
	 */
	BlockVector<double> prev_solution;
	/// The system right hand side vector.
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
		// material_id E nu alpha k
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
		// Command not found by this class
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

	// Set the boundaries pointer to that of the Executive class
	boundaries = bound;

	// Setup the system intially
	setup_system ();

	// Output the which time step is being analyzed
	if (is_transient)
		std::cout << "\n\nAnalyzing time: " << tc.get_current_time() << "\n";

	assemble_system ();
	solve ();
	output_results ();

	if (is_transient) {
		while (tc.increment_time()) {
			// Output the which time step is being analyzed
			std::cout << "\n\nAnalyzing time: " << tc.get_current_time() << "\n";
			// Call the transient assembly method, which doesn't update everything, only what is needed
			transient_assemble_system ();
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

	// Distribute the dof's to the FESystem
	dof_handler.distribute_dofs (fe);
	// Number the dof's component wise (be aware changing this renumbering will affect several parts of the code)
	DoFRenumbering::component_wise (dof_handler);

	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << "\n";

	// Find out how many dofs there are per component, hard coded that here is 4 components.  This is the case throughout
	// the code.  If the formulation changes, this will need to be looked at very carefully.
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

	// Set boundary id's
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
	}

	Status("Completed setup_system.", verbosity, MIN_V);
} // setup_system

// Private method: assemble_system
template<int dim>
void ThermalElasticityProblem<dim>::assemble_system()
{
	Status("Starting assemble_system.", verbosity, MIN_V);

	Status("Starting initialization of assembly variables.", verbosity, MAX_V);
	// Use Guass-Legendre quadrature of order 2 (2 quadrature points in each space direction), i.e. linear hexahedral elements
	QGauss<dim>  quadrature_formula(2);
	// By the nature of geometry, the face quadrature is a dimension lower than the problem
	QGauss<dim-1> face_quadrature_formula(2);

	// Extractors get the FESystem for a subset of the problem (used for shape functions mostly)
	const FEValuesExtractors::Scalar t_extract(0);
	const FEValuesExtractors::Vector u_extract(1);

	// Create FEValues and tell it what to update for points and faces
	FEValues<dim> fe_values (fe, quadrature_formula,
			update_values   | update_gradients |
			update_quadrature_points | update_JxW_values);
	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
			update_values         | update_quadrature_points  |
			update_normal_vectors | update_JxW_values);

	// Some variables that will be used a lot during the assembly
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();

	// Local matrix and vector of the cell. Though the problem's matrix is sparse, the local cell's matrix
	// is dense.
	FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>       cell_rhs (dofs_per_cell);

	// Vector used to map the local indicies to the global systems indicies
	std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	// Right hand side objects that contain heat generation and body force terms
	const ThermalRightHandSide<dim> thermal_rhs;
	const ElasticityRightHandSide<dim> elastic_rhs;
	std::vector<Vector<double> > elastic_rhs_values (n_q_points, Vector<double>(dim));

	Status("Completed initialization of assembly variables.", verbosity, MAX_V);

	Status("Starting the cell loop in assembly.", verbosity, MAX_V);

	// Iterate over every cell
	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		// Initialize the FEValues to the cell so the corresponding shape functions can be acquired
		fe_values.reinit (cell);

		// Make sure the cell's local matrix and rhs vector are zero since we only add terms (as opposed to setting terms)
		cell_matrix = 0;
		cell_rhs = 0;

		// Get global indicies corresponding to the cell's local indicies
		cell->get_dof_indices (local_dof_indices);

		// Copy the right hand side terms for the elasticity to a vector to be used in this cell's assembly
		elastic_rhs.vector_value_list (fe_values.get_quadrature_points(), elastic_rhs_values);

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
			// Loop over each dof of the cell (for each quad point)
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				// Temperature shape function for i
				const double phi_i_t = fe_values[t_extract].value(i, q_point);
				// Gradient of the temperature shape function for i
				const Tensor<1,dim>  grad_phi_i_t = fe_values[t_extract].gradient(i, q_point);
				// Displacement shape function for i
				const Tensor<1,dim>  phi_i_u = fe_values[u_extract].value(i, q_point);
				// Gradient of the displacement shape function for i (This is a symmetric tensor by nature and is stored
				// as a symmetric tensor since it makes operations with the material tensors easier)
				const SymmetricTensor<2,dim>  grad_phi_i_u = fe_values[u_extract].symmetric_gradient(i, q_point);
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					// Temperature shape function for j
					const double phi_j_t = fe_values[t_extract].value(j, q_point);
					// Gradient of the temperature shape function for j
					const Tensor<1,dim>  grad_phi_j_t = fe_values[t_extract].gradient(j, q_point);
					// Displacement shape function for j
					const SymmetricTensor<2,dim>  grad_phi_j_u = fe_values[u_extract].symmetric_gradient(j, q_point);

					/**
					 * The assembly of non-transient matrix terms:
					 * \f$(\nabla\phi_i^T*k*\nabla\phi_j^T - \nabla\phi_i^u*C*\alpha*\phi_j^T + \nabla\phi_i^u*C*\nabla\phi_j^u)*JxW\f$
					 */
					cell_matrix(i,j) += (grad_phi_i_t * k * grad_phi_j_t
							+ -1*(grad_phi_i_u * (stiffness * alpha) * phi_j_t)
							+ grad_phi_i_u * (stiffness * grad_phi_j_u)
							) *
							fe_values.JxW (q_point);
				} // j

				/**
				 * The assembly of non-transient rhs terms:
				 * \f$(\phi_i^T*Q + \phi_i^u*\rho*f)*JxW\f$
				 */
				cell_rhs(i) += (phi_i_t * thermal_rhs.value (fe_values.quadrature_point (q_point))
							+ phi_i_u * elastic_rhs_values[q_point]
							) * fe_values.JxW (q_point);
			} // i
		} // loop over quad points

		// Loop over faces and apply Nuemman BC's
		for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {
			// Check if face is at boundary first
			if (cell->face(face)->at_boundary()) {
				// Apply Nuemman BC for thermal flux
				for (unsigned int flux_bc_index = 0; flux_bc_index < flux_bcs.size(); flux_bc_index++) {
					if (cell->face(face)->boundary_indicator() == flux_bcs[flux_bc_index]->get_id()) {
						fe_face_values.reinit(cell, face);
						for (unsigned int q_point = 0; q_point<n_face_q_points; ++q_point) {
							double flux_value = flux_bcs[flux_bc_index]->value(fe_face_values.quadrature_point(q_point));

							// Modify the right hand side for the boundary condition
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

							// Modify the right hand side for the boundary condition
							for (unsigned int i=0; i<dofs_per_cell; ++i) {
								// Get the component of i
								const unsigned component_i = fe.system_to_component_index(i).first;
								// Only add the traction contribution if it is in the component mask
								if (traction_bcs[traction_bc_index]->get_comp_mask()[component_i] == true) {
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
			} // if cell at boundary
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
	Status("Completed the cell loop in assembly.", verbosity, MAX_V);

	// Apply hanging node constraints
	hanging_node_constraints.condense(system_matrix);
	hanging_node_constraints.condense(system_rhs);

	// Apply dirchlet BCs
	apply_dirchlet_bcs ();

	Status("Completed assemble_system.", verbosity, MIN_V);
} // assemble_system

// Private method: assemble_system
template<int dim>
void ThermalElasticityProblem<dim>::transient_assemble_system()
{
	// The documentaion of this method closely follows assemble_system(), see that method for clarification if needed

	// At each new time set, set the thermal solution and system block to 0 since terms would be continued to be added for
	// each time step
	system_matrix.block(0, 0) = 0;
	system_rhs.block(0) = 0;

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
	Status("Completed initialization of assembly variables.", verbosity, MAX_V);

	Status("Starting the cell loop in assembly.", verbosity, MAX_V);
	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{

		fe_values.reinit (cell);
		cell_matrix = 0;
		cell_rhs = 0;
		cell->get_dof_indices (local_dof_indices);

		// If transient get the previous solution
		if (is_transient) {
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				local_prev_solution[i] = prev_solution(local_dof_indices[i]);
		}

		// Get cell's material properties
		Tensor<2, dim> k;

		bool mat_found = false;
		for (unsigned int mat_ind = 0; mat_ind < materials.size(); mat_ind++)
			if (materials[mat_ind].get_id() == cell->material_id()) {
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
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					// Temperature shape function for j
					const double phi_j_t = fe_values[t_extract].value(j, q_point);
					const Tensor<1,dim>  grad_phi_j_t = fe_values[t_extract].gradient(j, q_point);

					// Only assemble terms corresponding to the thermal block
					cell_matrix(i,j) += (grad_phi_i_t * k * grad_phi_j_t) * fe_values.JxW (q_point);

					// Add transient terms
					cell_matrix(i,j) += (1/(tc.get_current_time() - tc.get_prev_time())) * phi_i_t * phi_j_t * fe_values.JxW (q_point);
					cell_rhs(i) += (1/(tc.get_current_time() - tc.get_prev_time())) * phi_i_t * phi_j_t * local_prev_solution[j] * fe_values.JxW (q_point);
				}

				// Only assembly right side terms corresponding to the thermal block
				cell_rhs(i) += (phi_i_t * thermal_rhs.value (fe_values.quadrature_point (q_point))) * fe_values.JxW (q_point);
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
	Status("Completed the cell loop in assembly.", verbosity, MAX_V);

	hanging_node_constraints.condense(system_matrix);
	hanging_node_constraints.condense(system_rhs);

	apply_dirchlet_bcs ();

	Status("Completed assemble_system.", verbosity, MIN_V);
} // assemble_system

template <int dim>
void ThermalElasticityProblem<dim>::apply_dirchlet_bcs ()
{
	// TODO: Consider deriving rhs from TensorFunction rather than just Function,
	// since it is more logical to describe as tensors

	// TODO: the component mask is stored in the VectorBoundary, allowing the boundary to only apply
	// to a specific component and not all components for a boundary id, but the input script has not
	// been modified to take in the values like that.  The code is ready for it, so if need just change
	// the input script reader.

	// Apply Dirchlet BC for temperature
	for (unsigned int i = 0; i < temperature_bcs.size(); i++) {
		std::map<unsigned int,double> boundary_values_map;
		VectorTools::interpolate_boundary_values (dof_handler,
					temperature_bcs[i]->get_id(),
					*temperature_bcs[i],
					boundary_values_map,
					temperature_bcs[i]->get_comp_mask());

		MatrixTools::apply_boundary_values (boundary_values_map,
			  system_matrix,
			  solution,
			  system_rhs);
	}

	// Apply Dirchlet BC for displacement
	for (unsigned int i = 0; i < displacement_bcs.size(); i++) {
		std::map<unsigned int,double> boundary_values_map;
		VectorTools::interpolate_boundary_values (dof_handler,
					displacement_bcs[i]->get_id(),
					*displacement_bcs[i],
					boundary_values_map,
					displacement_bcs[i]->get_comp_mask());

		MatrixTools::apply_boundary_values (boundary_values_map,
			  system_matrix,
			  solution,
			  system_rhs);
	}
}

template <int dim>
void ThermalElasticityProblem<dim>::solve ()
{
	Status("Starting solve for temperature.", verbosity, MIN_V);

	// We are just using the conjugate gradient method for this problem
	SolverControl           solver_control (1000, 1e-14);
	SolverCG<>              cg (solver_control);

	// If the problem is transient, the initial temperature is set and will be taken as the thermal solution at time step 0
	if (!is_transient || tc.index() > 0) {
		// First solve for the temperatures
		cg.solve(system_matrix.block(0, 0), solution.block(0), system_rhs.block(0), PreconditionIdentity());
		hanging_node_constraints.distribute(solution.block(0));
	}

	// Now solve for the displacements
	// Multiply the thermal solution with coupling block and subtract from rhs, but use the tmp Vector so the
	// actual rhs vector is not modified, allowing it to stay the same across time steps.
	Vector<double> tmp (solution.block(1).size());
	system_matrix.block(1, 0).vmult(tmp, solution.block(0));
	tmp *= -1;
	tmp += system_rhs.block(1);

	Status("Starting solve for displacement.", verbosity, MIN_V);

	PreconditionSSOR<> preconditioner;
	preconditioner.initialize(system_matrix.block(1, 1), 1.2);
	cg.solve(system_matrix.block(1, 1), solution.block(1), tmp, preconditioner);
	hanging_node_constraints.distribute(solution.block(1));

	//Status("Completed solve for displacement.", verbosity, MIN_V);

	if (is_transient) {
		Status("Saving current solution for next time step.", verbosity, MIN_V);
		prev_solution = solution;
	}
} // solve

template <int dim>
void ThermalElasticityProblem<dim>::output_results () const
{

	// Ok, so here things get a bit complicated since the PostProcesor class does not allow the compute_derived_quantities
	// method to know which cell it is operatoring on or the material id, which is necessary to compute the stress.
	// See the TEPostProcessor and TEDataOut methods for more info on what is going on in there, but here just treat
	// it like a normal PostProcessor and DataOut class
	TEPostProcessor<dim> post_process(materials);
	TEDataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, post_process);
	// Note to self here: an arguement to this method will subdived the cells and interpolate the values.  This is very
	// useful if quadratic elements are used, but since we use linear here this is less interesting.
	data_out.build_patches();
	// Pick the name depending on if the problem is transient or not
	// TODO: The input script should have a parameter for the output file name
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
