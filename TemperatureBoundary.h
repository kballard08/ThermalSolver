/*
 * TemperatureBoundary.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef TEMPERATUREBOUNDARY_H_
#define TEMPERATUREBOUNDARY_H_

#include "BoundaryCondition.h"

template<int dim>
class TemperatureBoundary : public BoundaryCondition
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


#endif /* TEMPERATUREBOUNDARY_H_ */
