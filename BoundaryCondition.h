/*
 * BoundaryCondition.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <deal.II/base/point.h>

using namespace dealii;

// The BoundaryCondition class is meant to describe a boundary condition for the problem.  It is an abstract class and the inherited classes
// implement the class for different types of BC's.

// Use the PointInBoudnary method to check if a point is within the geometrical description
// of the boudnary's domain.  The domain of the boudnary is very vague is very flexible.  It just needs to cover the domain of the applicable
// boundary without touching any other boundaries, but it is fine if the domain covers not boundary areas.  For example, if the boundary applies
// to all boundary points within the cube (0, 0, 0) to (1, 1, 1), the PointInBoundary will evaluate the argument to see if it falls within that
// range and return the boolean result.

// The Value method is used to get the actual value of the boundary value for a given point.

template<int dim>
class BoundaryCondition
{
public:
	BoundaryCondition() {};
	virtual ~BoundaryCondition() {};

	// Virtual Methods
	virtual bool	point_in_boundary(const Point<dim> &p) = 0;
	virtual double	value(const Point<dim> &p) = 0;
};


#endif /* BOUNDARYCONDITION_H_ */
