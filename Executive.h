/*
 * Executive.h
 *
 *  Created on: Mar 7, 2013
 *      Author: kballard
 */

#ifndef EXECUTIVE_H_
#define EXECUTIVE_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "ElasticityProblem.h"
#include "ThermalProblem.h"
#include "BoundaryGeometry.h"
#include "ScriptReader.h"

namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class Executive
{
public:
	Executive();
	~Executive();

	void run(ScriptReader *script_reader);

private:
	void make_grid_test ();

	Verbosity verbosity;

	Triangulation<dim>  triangulation;
	std::vector<BoundaryGeometry<dim> *> *boundaries;
	ScriptReader* 		sr;

	ElasticityProblem<dim>	elasticity_problem;
	ThermalProblem<dim>		thermal_problem;

	// State variables
	bool mesh_initialized;
};

// Constructor
template<int dim>
Executive<dim>::Executive() : elasticity_problem(&triangulation), thermal_problem(&triangulation)
{
	// Make sure the analysis only exists for 2d and 3d
	Assert (dim == 2 || dim == 3, ExcNotImplemented());

	verbosity = MIN_V;
	boundaries = new std::vector<BoundaryGeometry<dim> *>();
	sr = 0;

	mesh_initialized = false;
}

// Destructor
template<int dim>
Executive<dim>::~Executive()
{
	// If the pointer is valid clean up
	if (boundaries != 0) {
		for (unsigned int i = 0; i < boundaries->size(); i++) {
			delete (*boundaries)[i];
		}
		delete boundaries;
	}
}

// Public member: run
template <int dim>
void Executive<dim>::run(ScriptReader *script_reader)
{
	sr = script_reader;

	std::vector<std::string> tokens;
	int lineNum = 0;
	while(sr->get_next_line(tokens)) {
		if (verbosity >= MAX_V) {
			std::cout << "Processing line: " << lineNum++;
			for (unsigned int i = 0; i < tokens.size(); i++)
				std::cout << " " <<tokens[i];
			std::cout << std::endl;
		}

		// TODO: make a more efficient process command structure
		// Process the commands of the input file
		if (tokens[0] == "UseDefaultMesh") {
			make_grid_test();
		} // TODO: add other mesh/grid methods to read from file
		else if (tokens[0] == "SetVerbosity") {
			Assert(tokens.size() == 2, ExcMessage("The input line SetVerbosity expects an argument for the verbosity level to set."))
			verbosity = (Verbosity)atoi(tokens[1].c_str());
		}
		else if (tokens[0] == "ReadBoundaries") { // ReadBoundaries
			while (sr->get_next_line(tokens)) {
				if (verbosity >= MAX_V) {
					std::cout << "Processing line: " << lineNum++;
					for (unsigned int i = 0; i < tokens.size(); i++)
						std::cout << " " <<tokens[i];
					std::cout << std::endl;
				}

				if (tokens[0] == "EndBoundaries") {
					Status("Finished reading boundaries.", verbosity, MAX_V);
					break;
				}
				else if (tokens.size() == dim*dim + 1) {
					// For 2d: (forms a line and assigns boundary id to line)
					// boundary_id x1 y1 x2 y2
					// For 3d: (forms a line and assigns boundary id to line)
					// boundary_id x1 y1 z1 x2 y2 z2 x3 y3 z3
					BoundaryGeometry<dim> * bg_p = new BoundaryGeometry<dim>(atoi(tokens[0].c_str()));

					// Get the points
					int index = 1;
					Point<dim> points[dim];
					for (int i = 0; i < dim; i++)
						for (int j = 0; j < dim; j++) {
							points[i](j) = atof(tokens[index].c_str());
							index++;
						}

					bg_p->define_geometry(points);
					boundaries->push_back(bg_p);
				}
				else if (tokens.size() == 3) {
					// Expected tokens for 2d and 3d are:
					// boundary_id axis value_along_axis
					BoundaryGeometry<dim> * bg_p = new BoundaryGeometry<dim>(atoi(tokens[0].c_str()));
					bg_p->define_geometry(atoi(tokens[1].c_str()), atof(tokens[2].c_str()));
					boundaries->push_back(bg_p);
				}
				else {
					Assert(false, ExcNotImplemented())
				}
			}
		} // ReadBoundaries
		else if (tokens[0] == "ReadBCs") { // ReadBCs
			while (sr->get_next_line(tokens)) {
				if (tokens[0] == "EndBCs") {
					Status("Finished reading boundary conditions.", verbosity, MAX_V);
					break;
				}
				else {
					bool process_result = false;
					process_result = elasticity_problem.process_bc(tokens);
					if (process_result == false)
						process_result = thermal_problem.process_bc(tokens);
					if (process_result == false)
						Assert(false, ExcNotImplemented())
				}
			}
		} // ReadBCs
		else { // If command is not recognized by executive, pass to elasticity and thermal problems to see if they can use it
			bool process_result = false;
			process_result = elasticity_problem.process_cmd(tokens);
			if (process_result == false)
				process_result = thermal_problem.process_cmd(tokens);
			if (process_result == false)
				Assert(false, ExcNotImplemented()) // The executive, elasticity problem, and thermal problem cannot use it
		}
	}

	// Let the thermal problem solve
	thermal_problem.run(boundaries);

	// Solve the elasticity problem
	elasticity_problem.run(boundaries);
}

// Private method: make_grid
template<int dim>
void Executive<dim>::make_grid_test()
{
	Status("Using the test grid.", verbosity, MIN_V);

	// For now just generate cube
	// Later include functionality to read in a mesh file?
	GridGenerator::hyper_cube(triangulation, -1, 1);
	triangulation.refine_global(4);

	// Update state information
	mesh_initialized = true;
}

}

#endif /* EXECUTIVE_H_ */
