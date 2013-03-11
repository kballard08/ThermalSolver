/*
 * BoundaryCondition.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>

namespace FEASolverNS
{

using namespace dealii;

// The BoundaryCondition class is meant to describe a boundary condition for the problem.  It is an abstract class and the inherited classes
// implement the class for different types of BC's.

// The Value method is used to get the actual value of the boundary value for a given point.

template<int dim>
class BoundaryCondition : public Function<dim>
{
public:
	BoundaryCondition() : Function<dim>() {};
	virtual ~BoundaryCondition() {};

	// Virtual Methods
	virtual int get_id() = 0;
	virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const = 0;
};

}

#endif /* BOUNDARYCONDITION_H_ */
