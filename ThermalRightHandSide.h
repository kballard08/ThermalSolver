/*
 * RightHandSide.h
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#ifndef RIGHTHANDSIDE_H_
#define RIGHTHANDSIDE_H_

#include <deal.II/base/function.h>

namespace FEASolverNS
{

using namespace dealii;

template <int dim>
class ThermalRightHandSide : public Function<dim>
{
public:
	ThermalRightHandSide() {};
	virtual ~ThermalRightHandSide() {};

	virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
inline
double ThermalRightHandSide<dim>::value (const Point<dim> &p, const unsigned int /*component*/) const
{
	// For now this assumes that there is no heat generation in the body
	// which is fine for now.  Later this may be a function of position
	// and time
	return 0;
}

}

#endif /* RIGHTHANDSIDE_H_ */
