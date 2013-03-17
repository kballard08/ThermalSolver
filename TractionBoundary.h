/*
 * TractionBoundary.h
 *
 *  Created on: Mar 13, 2013
 *      Author: kballard
 */

#ifndef TRACTIONBOUNDARY_H_
#define TRACTIONBOUNDARY_H_

#include "BoundaryCondition.h"
#include <vector>

namespace FEASolverNS
{

template<int dim>
class TractionBoundary : public BoundaryCondition<dim>
{
public:
	// Constructor
	TractionBoundary(int boundary_id, const std::vector<double> &traction);
	virtual ~TractionBoundary() {};

	virtual int get_id() { return bound_id; };

	virtual double value(const Point<dim> &p, const unsigned int component = 0) const;

	virtual void vector_value(const Point<dim> &p, Vector<double>   &values) const;

	virtual void vector_value_list(const std::vector<Point<dim>> &points, std::vector<Vector<double>>   &value_list) const;

private:
	int bound_id;
	std::vector<double> tract;
};

template <int dim>
TractionBoundary<dim>::TractionBoundary(int boundary_id, const std::vector<double> &traction) : bound_id(boundary_id), tract(traction)
{
	// Make sure the prescribed traction is of the correct dimension
	Assert(traction.size() == dim, ExcDimensionMismatch (traction.size(), dim));
}

template <int dim>
double TractionBoundary<dim>::value(const Point<dim> &p, const unsigned int component) const
{
	return tract[component];
}

template <int dim>
void TractionBoundary<dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
	// Check for dimension mismatch between the force vector (values) and the dimension
	Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));

	// Just deep copy the tract vector
	for (int i = 0; i < dim; i++)
		values(i) = tract[i];
}

template <int dim>
void TractionBoundary<dim>::vector_value_list (const std::vector<Point<dim>> &points, std::vector<Vector<double>> &value_list) const
{
	// Make sure the number of points match the number of value vectors passed
	Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));

	// Number of points
	const unsigned int n_points = points.size();

	// Use the vector_value method on each point
    for (unsigned int p=0; p<n_points; ++p)
    	TractionBoundary<dim>::vector_value(points[p], value_list[p]);
}

}

#endif /* TRACTIONBOUNDARY_H_ */
