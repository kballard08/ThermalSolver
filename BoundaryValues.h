/*
 * BoundaryValues.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef BOUNDARYVALUES_H_
#define BOUNDARYVALUES_H_

#include <deal.II/base/function.h>
#include <deal.II/base/exceptions.h>
#include <vector>
#include <iostream>

#include "BoundaryCondition.h"

namespace ThermalSolverNS
{

using namespace dealii;

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
	BoundaryValues() {};
	virtual ~BoundaryValues();

	void add_boundary_condition(BoundaryCondition<dim> * bc_p) { boundary_conditions.push_back(bc_p); };

	virtual double value (const Point<dim>   &p,
						const unsigned int  component = 0) const;

private:
	// Stores pointers to the boundary conditions
	std::vector< BoundaryCondition<dim> * > boundary_conditions;
};

// Destructor
template <int dim>
inline
BoundaryValues<dim>::~BoundaryValues()
{
	// Clean up the vector of pointers
	for(unsigned int i = 0; i < boundary_conditions.size(); i++)
		delete boundary_conditions[i];
}

// Method: value
template <int dim>
inline
double BoundaryValues<dim>::value (const Point<dim> &p, const unsigned int /*component*/) const
{
	// Ouput the point for debugging
	//std::cout << "Finding Point (";
	//for (int j = 0; j < dim; j++)
	//	std::cout << p(j) << " ";
	//std::cout << ")" << std::endl;

	// Loop through the boundary conditions and fine the one that contains the point
	for(unsigned int i = 0; i < boundary_conditions.size(); i++)
	{
		if (boundary_conditions[i]->point_in_boundary(p))
			return boundary_conditions[i]->value(p);
	}

	Assert(false, ExcMessage("A boundary point was not found in the boundary conditions."))

	// Return 0 so the compiler will be happy that every path has a return value, but this will never be reached.
	return 0;
}

}

#endif /* BOUNDARYVALUES_H_ */
