/*
 * FluxBoundary.h
 *
 *  Created on: Mar 4, 2013
 *      Author: kballard
 */

#ifndef FLUXBOUNDARY_H_
#define FLUXBOUNDARY_H_

#include "BoundaryCondition.h"

namespace ThermalSolverNS
{

template<int dim>
class FluxBoundary : public BoundaryCondition<dim>
{
public:
	// Constructor
	// Forms a rectangular prism using the two points for the domain of the boundary
	FluxBoundary(int boundary_id, double val) : bound_id(boundary_id), flux_value(val) {};
	virtual ~FluxBoundary() {};

	virtual int get_id() { return bound_id; };
	virtual double value (const Point<dim> &p, const unsigned int component = 0) const;

private:
	int bound_id;
	double flux_value;
};

template <int dim>
double FluxBoundary<dim>::value(const Point<dim> &p, const unsigned int component) const
{
	// For now assume it is a static temperature
	// In the future this could vary with location
	// or time (in transient problems)
	return flux_value;
}

}

#endif /* FLUXBOUNDARY_H_ */
