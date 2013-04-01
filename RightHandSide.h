/*
 * RightHandSide.h
 *
 *  Created on: Mar 13, 2013
 *      Author: kballard
 */

#ifndef RIGHTHANDSIDE_H_
#define RIGHTHANDSIDE_H_

#include <deal.II/base/function.h>

#include <vector>

namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class RightHandSide :  public Function<dim>
{
public:
	RightHandSide () : Function<dim> (dim+1) {};
	virtual ~RightHandSide() {};

	virtual void vector_value (const Point<dim> &p,
				Vector<double>   &values) const;

	virtual void vector_value_list (const std::vector<Point<dim> > &points,
				std::vector<Vector<double> >   &value_list) const;
};

template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
	// Check for dimension mismatch between the force vector (values) and the dimension
	Assert (values.size() == dim+1, ExcDimensionMismatch (values.size(), dim+1));
	// Make sure the analysis only exists for 2d and 3d
	Assert (dim == 2 || dim == 3, ExcNotImplemented());

	// For now, assume no body forces exist and no heat generation
	// DoF 0 is for temperature
	// DoF 1-dim is for displacement
	// Later this can be a function of position or material or something
	for (int i = 0; i < dim+1; i++)
		values(i) = 0;
}

template <int dim>
inline
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim>> &points, std::vector<Vector<double>> &value_list) const
{
	// Make sure the number of points match the number of value vectors passed
	Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));

	// Number of points
	const unsigned int n_points = points.size();

	// Use the vector_value method on each point
    for (unsigned int p=0; p<n_points; ++p)
    	RightHandSide<dim>::vector_value(points[p], value_list[p]);
}

}


#endif /* RIGHTHANDSIDE_H_ */
