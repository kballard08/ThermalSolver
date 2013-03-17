/*
 * DisplacementBoundary.h
 *
 *  Created on: Mar 11, 2013
 *      Author: kballard
 */

#ifndef DISPLACEMENTBOUNDARY_H_
#define DISPLACEMENTBOUNDARY_H_

#include "BoundaryCondition.h"
#include <vector>

namespace FEASolverNS
{

template<int dim>
class DisplacementBoundary : public BoundaryCondition<dim>
{
public:
	// Constructor
	DisplacementBoundary(int boundary_id, const std::vector<double> &displacement);
	virtual ~DisplacementBoundary() {};

	virtual int get_id() { return bound_id; };

	virtual double value (const Point<dim> &p, const unsigned int component = 0) const;

	virtual void vector_value (const Point<dim> &p, Vector<double>   &values) const;

	virtual void vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> >   &value_list) const;

private:
	int bound_id;
	std::vector<double> disp;
};

template <int dim>
DisplacementBoundary<dim>::DisplacementBoundary(int boundary_id, const std::vector<double> &displacement) : bound_id(boundary_id), disp(displacement)
{
	// Make sure the prescribed displacement is of the correct dimension
	Assert(displacement.size() == dim, ExcDimensionMismatch (displacement.size(), dim));
}

template <int dim>
double DisplacementBoundary<dim>::value(const Point<dim> &p, const unsigned int component) const
{
	return disp[component];
}

template <int dim>
void DisplacementBoundary<dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
	// Check for dimension mismatch between the force vector (values) and the dimension
	Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));

	// Just deep copy the disp vector
	for (int i = 0; i < dim; i++)
			values(i) = disp[i];
}

template <int dim>
void DisplacementBoundary<dim>::vector_value_list (const std::vector<Point<dim>> &points, std::vector<Vector<double>> &value_list) const
{
	// Make sure the number of points match the number of value vectors passed
	Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));

	// Number of points
	const unsigned int n_points = points.size();

	// Use the vector_value method on each point
    for (unsigned int p=0; p<n_points; ++p)
    	DisplacementBoundary<dim>::vector_value(points[p], value_list[p]);
}

}


#endif /* DISPLACEMENTBOUNDARY_H_ */
