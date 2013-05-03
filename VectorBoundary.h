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
#include "deal.II/fe/component_mask.h"
#include <map>
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

/**
 * The VectorBoundary class describes a boundary condition of either Dirchlet or Nuemman types, and can be used for others I'm sure.
 * The idea is to pass the class the values and components that it applies to and the condition will be applied during assembly.
 * This is simply a storage class for all intents and purposes.  The vector that it is added to in the ThermalElsaticityProblem
 * class determines what type of boundary condition it is.
 */
class VectorBoundary : public Function<dim>
{
public:
	/// The constructor takes the boundary id (must be unique otherwise the last one added will be applied), the values, the component
	/// indicies it applies to, and the number of components of the FESystem
	VectorBoundary(const int &boundary_id, const std::vector<double> &values, const std::vector<unsigned int> &indices, const unsigned int &n_components);
	/// Virtual destructor
	virtual ~VectorBoundary() {};

	/// Gets the corresponding boundary id
	int get_id() { return bound_id; };
	/// Get the component mask for the components that were passed to the constructor
	ComponentMask get_comp_mask() { return ComponentMask(comp_mask); };
	/// Returns a map of components to value
	std::map<unsigned int, double> get_val_map() { return val_map; };
	/// Returns the vector of values from the start index to end index
	void vector_value (const Point<dim> &p, Vector<double>   &values, const unsigned int &start_ind, const unsigned int &end_ind) const;

	// Virtual Methods from function class
	virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
	virtual void vector_value (const Point<dim> &p, Vector<double>   &values) const;
	virtual void vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> >   &value_list) const;

	unsigned int val_size() const { return val.size(); };
	double val_value(unsigned int i) const { Assert(i >= 0 && i < val.size(), ExcIndexRange(i, 0, val.size())); return val[i]; };

protected:
	int bound_id;
	std::vector<double> val;
	std::map<unsigned int, double> val_map;
	std::vector<bool> comp_mask; // Internally store the component mask for use in applying the boundary conditionsss

};

template <int dim>
VectorBoundary<dim>::VectorBoundary(const int &boundary_id,
									const std::vector<double> &values,
									const std::vector<unsigned int> &indices,
									const unsigned int &n_components)
: Function<dim>(n_components), bound_id(boundary_id), val(n_components), comp_mask(n_components, false)
{
	// Make sure the prescribed values are of the correct dimension
	Assert(values.size() <= n_components, ExcMessage("The size of the value vector passed to the VectorBoundary must be <= n_components"));
	Assert(values.size() == indices.size(), ExcDimensionMismatch(values.size(), indices.size()));

	for (unsigned int i = 0; i < values.size(); i++) {
		// Make sure index is valid
		Assert(indices[i] >= 0 && indices[i] < n_components, ExcIndexRange(indices[i], 0, n_components));

		val[ indices[i] ] = values[i];
		comp_mask[ indices[i] ] = true;
		val_map[ indices[i] ] = values[i];
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
	for (unsigned int i = 0; i < val.size(); i++)
			values(i) = val[i];
}

template <int dim>
void VectorBoundary<dim>::vector_value (const Point<dim> &p, Vector<double>   &values, const unsigned int &start_ind, const unsigned int &end_ind) const
{
	unsigned int tmp = end_ind - start_ind + 1;

	// Check for dimensions
	Assert (values.size() == tmp, ExcDimensionMismatch (values.size(), tmp));
	Assert (start_ind >= 0 && start_ind <= end_ind, ExcIndexRange (start_ind, 0, end_ind+1));
	Assert (end_ind >= start_ind && end_ind < val.size(), ExcIndexRange (end_ind, start_ind, val.size()));

	// Copy indices [start_ind, end_ind] from val into [0, start_ind-end_ind+1] of values
	unsigned int values_ind = 0;
	for (unsigned int i = start_ind; i <= end_ind; i++, values_ind++)
			values(values_ind) = val[i];
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
