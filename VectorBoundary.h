/*
 * VectorBoundary.h
 *
 *  Created on: Mar 13, 2013
 *      Author: kballard
 */

#ifndef VECTOROUNDARY_H_
#define VECTORBOUNDARY_H_

#include "deal.II/base/function.h"
#include "deal.II/base/tensor_function.h"
#include <iostream>
#include <vector>

namespace FEASolverNS
{

using namespace dealii;

// TODO: Document this class
// Has dim+1 components but only uses dim values of it
// start_index indicates the beginning of the data in the dim+1
// Indices indicate the component indices in the dof_handler that the values correspond to

template<int dim>


class VectorBoundary : public Function<dim>
{
public:
	// Constructor
	VectorBoundary(const int &boundary_id, const std::vector<double> &values, const std::vector<unsigned int> &indices, const unsigned int &n_components);
	virtual ~VectorBoundary() {};

	virtual int get_id() { return bound_id; };

	virtual double value (const Point<dim> &p, const unsigned int component = 0) const;

	virtual void vector_value (const Point<dim> &p, Vector<double>   &values) const;

	virtual void vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> >   &value_list) const;

	virtual void nonzero_vector_value (const Point<dim> &p, Vector<double>   &values) const;

	unsigned int val_size() const { return val.size(); };
	double val_value(unsigned int i) const { Assert(i >= 0 && i < val.size(), ExcIndexRange(i, 0, val.size())); return val[i]; };

	//friend std::ostream& operator<<(std::ostream& os, const VectorBoundary<dim>& vb);
protected:
	int bound_id;
	std::vector<double> val;
	std::vector<double> nz_val; // nonzero values, the unmodified values vector passed to the constructor
	std::vector<bool> comp_mask; // Internally store the component mask for use in applying the boundary conditionsss


};

template <int dim>
VectorBoundary<dim>::VectorBoundary(const int &boundary_id,
									const std::vector<double> &values,
									const std::vector<unsigned int> &indices,
									const unsigned int &n_components)
: Function<dim>(n_components), bound_id(boundary_id), val(n_components), nz_val(values), comp_mask(n_components)
{
	// Make sure the prescribed values are of the correct dimension
	Assert(values.size() <= n_components, ExcMessage("The size of the value vector passed to the VectorBoundary must be <= n_components"));
	Assert(values.size() == indices.size(), ExcDimensionMismatch(values.size(), indices.size()));

	for (unsigned int i = 0; i < values.size(); i++) {
		// Make sure index is valid
		Assert(indices[i] >= 0 && indices[i] < n_components, ExcIndexRange(indices[i], 0, n_components));

		val[ indices[i] ] = values[i];
	}
}

template <int dim>
double VectorBoundary<dim>::value(const Point<dim> &p, const unsigned int component) const
{
	Assert(component >= 0 && component < val.size(), ExcIndexRange(component, 0, val.size()));
	return val[component];
}

template <int dim>
void VectorBoundary<dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
	// Check for dimension mismatch between the force vector (values) and the dimension
	Assert (values.size() == val.size(), ExcDimensionMismatch (values.size(), val.size()));

	// Just deep copy the val vector
	for (int i = 0; i < dim; i++)
			values(i) = val[i];
}

template <int dim>
void VectorBoundary<dim>::vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> > &value_list) const
{
	// Make sure the number of points match the number of value vectors passed
	Assert (value_list.size() == points.size(), ExcDimensionMismatch (value_list.size(), points.size()));

	// Number of points
	const unsigned int n_points = points.size();

	// Use the vector_value method on each point
    for (unsigned int p=0; p<n_points; ++p)
    	VectorBoundary<dim>::vector_value(points[p], value_list[p]);
}

template <int dim>
void VectorBoundary<dim>::nonzero_vector_value (const Point<dim> &p, Vector<double> &values) const
{
	// Check for dimension mismatch between the force vector (values) and the dimension
	Assert (values.size() == nz_val.size(), ExcDimensionMismatch (values.size(), nz_val.size()));

	// Just deep copy the val vector
	for (unsigned int i = 0; i < nz_val.size(); i++)
			values(i) = nz_val[i];
}

// Friend class for ostream output
template <int dim>
std::ostream& operator<<(std::ostream& os, const VectorBoundary<dim>& vb) {
	os << "[";
	for (unsigned int i = 0; i < vb.val_size(); i++) {
		os << vb.val_value(i);
		if (i != vb.val_size() - 1)
			os << ", ";
	}
	os << "]";
	return os;
}

}

#endif /* VECTOROUNDARY_H_ */
