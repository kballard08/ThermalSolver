/*
 * TemperatureBoundary.cpp
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#include "TemperatureBoundary.h"

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
