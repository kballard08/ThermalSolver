/*
 * ElasticityRightHandSide.h
 *
 *  Created on: Mar 13, 2013
 *      Author: kballard
 */

#ifndef ELASTICITYRIGHTHANDSIDE_H_
#define ELASTICITYRIGHTHANDSIDE_H_

#include <deal.II/base/function.h>

#include <vector>

namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class ElasticityRightHandSide :  public Function<dim>
{
public:
	ElasticityRightHandSide () : Function<dim> (dim) {};
	virtual ~ElasticityRightHandSide() {};

	virtual void vector_value (const Point<dim> &p,
				Vector<double>   &values) const;

	virtual void vector_value_list (const std::vector<Point<dim> > &points,
				std::vector<Vector<double> >   &value_list) const;
};

template <int dim>
inline
void ElasticityRightHandSide<dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
	// Check for dimension mismatch between the force vector (values) and the dimension
	Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));
	// Make sure the analysis only exists for 2d adn 3d
	Assert (dim == 2 || dim == 3, ExcNotImplemented());

	// For now, assume no body forces exist
	for (int i = 0; i < dim; i++)
		values(i) = 0;
}

template <int dim>
inline
void ElasticityRightHandSide<dim>::vector_value_list (const std::vector<Point<dim>> &points, std::vector<Vector<double>> &value_list) const
{
	// Make sure the number of points match the number of value vectors passed
	Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));

	// Number of points
	const unsigned int n_points = points.size();

	// Use the vector_value method on each point
    for (unsigned int p=0; p<n_points; ++p)
    	ElasticityRightHandSide<dim>::vector_value(points[p], value_list[p]);
}

}


#endif /* ELASTICITYRIGHTHANDSIDE_H_ */
