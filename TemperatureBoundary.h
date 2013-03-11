/*
 * TemperatureBoundary.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef TEMPERATUREBOUNDARY_H_
#define TEMPERATUREBOUNDARY_H_

#include "BoundaryCondition.h"

namespace FEASolverNS
{

template<int dim>
class TemperatureBoundary : public BoundaryCondition<dim>
{
public:
	// Constructor
	// Forms a rectangular prism using the two points for the domain of the boundary
	TemperatureBoundary(int boundary_id, double T0) : bound_id(boundary_id), Temp0(T0) {};
	virtual ~TemperatureBoundary() {};

	virtual int get_id() { return bound_id; };
	virtual double value (const Point<dim> &p, const unsigned int component = 0) const;

private:
	int bound_id;
	double Temp0;
};

template <int dim>
double TemperatureBoundary<dim>::value(const Point<dim> &p, const unsigned int component) const
{
	// For now assume it is a static temperature
	// In the future this could vary with location
	// or time (in transient problems)
	return Temp0;
}

}

#endif /* TEMPERATUREBOUNDARY_H_ */
