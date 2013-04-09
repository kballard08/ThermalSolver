/*
 * BoundaryGeometry.h
 *
 *  Created on: Feb 27, 2013
 *      Author: kballard
 */

#ifndef BOUNDARYGEOMETRY_H_
#define BOUNDARYGEOMETRY_H_

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <math.h>

namespace FEASolverNS
{

using namespace dealii;

// Plane for 3d and line for 2d
template<int dim>
class BoundaryGeometry
{
public:
	BoundaryGeometry(int boundary_id) : bound_id(boundary_id), defined_by_tensor(false), defined_by_axis(false), def_axis(0), def_value(0) {};
	~BoundaryGeometry() {};

	// Pass the points to define the geometry
	// 2d: 2 points (line)
	// 3d: 3 points (plane)
	void define_geometry(Point<dim> points[dim]);
	// Defines a geometry normal to a principle axis at a given value along the axis
	// 0: x
	// 1: y
	// 2: z
	void define_geometry(int axis, double value);

	// Test if point lies on the plane
	bool point_on_geometry(const Point<dim> &p);

	// Returns the associated boundary_id
	int get_id() { return bound_id; }

private:
	int bound_id;

	// Options to define
	bool defined_by_tensor;
	bool defined_by_axis;

	// Tensor representing the points forming the geometry
	Tensor<2, dim> def_points;

	int def_axis;
	double def_value;

	// Tweak this value if points that should be on a plane are not or vice versa
	const static double tol = 1e-12;
};

template<int dim>
void BoundaryGeometry<dim>::define_geometry(Point<dim> points[dim]) {
	// Tensore formed by:
	// [p1(x) p1(y) p1(z);
	//  p2(x) p2(y) p2(z);
	//  p3(x) p3(y) p3(z)]
	for (int i = 0; i < dim; i ++)
		for (int j = 0; j < dim; j++)
			def_points[i][j] = points[i](j);

	defined_by_tensor = true;
}

template<int dim>
void BoundaryGeometry<dim>::define_geometry(int axis, double value) {
	Assert(axis < dim && axis >= 0, ExcMessage("The axis passed to define_geometry must be within the analysis's dimensions."))

	def_axis = axis;
	def_value = value;
	defined_by_axis = true;
}

template<int dim>
bool BoundaryGeometry<dim>::point_on_geometry(const Point<dim> &p)
{
	if (defined_by_tensor) {
		// Make copy of tensor, so when we do some math we don't mess it up
		Tensor<2, dim> tmp_tensor = def_points;

		// Loop through and modify tensor to be:
		// [p(x)-p1(x) p(y)-p1(y) p(z)-p1(z);
		//  p(x)-p2(x) p(y)-p2(y) p(z)-p2(z);
		//  p(x)-p3(x) p(y)-p3(y) p(z)-p3(z)]
		for (int i = 0; i < dim; i ++)
			for (int j = 0; j < dim; j++)
				tmp_tensor[i][j] -= p(j);

		// If det(tmp_tensor) is within the set tolerance, the point is
		// understood to be on the plane (numerically of course)
		if (fabs(determinant(tmp_tensor)) <= tol)
			return true;
		else
			return false;
	}
	else if (defined_by_axis) {
		if (fabs(def_value - p(def_axis)) <= tol)
			return true;
		else
			return false;
	}
	else {
		Assert(false, ExcMessage("BoundaryGeometry must have define_geometry called before point_on_geometry is called."))
		return false; // Simply to make the compiler happy about return paths, it will never reach this point
	}
}

}

#endif /* BOUNDARYGEOMETRY_H_ */
