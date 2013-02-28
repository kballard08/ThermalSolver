/*
 * TemperatureBoundary.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef TEMPERATUREBOUNDARY_H_
#define TEMPERATUREBOUNDARY_H_

#include "BoundaryCondition.h"

namespace ThermalSolverNS
{

template<int dim>
class TemperatureBoundary : public BoundaryCondition<dim>
{
public:
	// Constructor
	// Forms a rectangular prism using the two points for the domain of the boundary
	TemperatureBoundary(Point<dim> p1, Point<dim> p2, double T0) : point1(p1), point2(p2), Temp0(T0) {};
	~TemperatureBoundary() {};

	bool	point_in_boundary(const Point<dim> &p);
	double	value(const Point<dim> &p);

private:
	Point<dim> point1, point2;
	double Temp0;
};

template <int dim>
bool TemperatureBoundary<dim>::point_in_boundary(const Point<dim> &p)
{
	// Loop over each component of the point (determined by the dimension)
	for (int i = 0; i < dim; i++)
	{
		if (point1(i) > point2(i)) // If coordinate i of the point1 is > point2
		{
			if (p(i) > point1(i) || p(i) < point2(i))
				return false; // p is outside bounds
		}
		else // If coordinate i of point 1 is <= point2
		{
			if (p(i) < point1(i) || p(i) > point2(i))
				return false; // p is outside bounds
		}
	}
	return true;
}

template <int dim>
double TemperatureBoundary<dim>::value(const Point<dim> &p)
{
	// For now assume it is a static temperature
	// In the future this could vary with location
	// or time (in transient problems)
	return Temp0;
}

}

#endif /* TEMPERATUREBOUNDARY_H_ */
