/*
 * TMPostProcessor.h
 *
 *  Created on: Apr 22, 2013
 *      Author: kballard
 */

#ifndef TEPOSTPROCESSOR_H_
#define TEPOSTPROCESSOR_H_

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/numerics/data_postprocessor.h>

#include "Material.h"

#include <string>
#include <vector>

namespace FEASolverNS
{

using namespace dealii;

/**
 * The TEPostProcessor class is intended to be paired with the TEDataOut class.
 * If it is used with the DataOut class, an exception will be thrown. It must
 * be specific like this to accomdate the output of stress data.  The major
 * change is that the compute_derived_quanitities_vector now accepts the
 * cell's material_id as an argument.
 */
template <int dim>
class TEPostProcessor : public DataPostprocessor<dim>
{
public:
	/// Constructor takes a reference to the materials vector, which it copies.  This should be made a pointer later, but
	/// since the memory for materials is very small for my cases, I am less worried about the performance hit.
	TEPostProcessor(const std::vector< Material<dim> > &materials) : mats(materials) {};

	/// The expected overload of compute_derived_quantities_vector from the DataPostprocessor class.
	/// Right now, the method will result in an Assertion explaining that it should never be called and
	/// to use the TEDataOut class with this class instead.
	virtual void compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
												   const std::vector<std::vector<Tensor<1,dim> > > &duh,
												   const std::vector<std::vector<Tensor<2,dim> > > &dduh,
												   const std::vector<Point<dim> >                  &normals,
												   const std::vector<Point<dim> >                  &evaluation_points,
												   std::vector<Vector<double> >                    &computed_quantities) const;

	/// Modified compute_derived_quantities_vector that takes the material id of the cell as well
	void compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
												   const std::vector<std::vector<Tensor<1,dim> > > &duh,
												   const std::vector<std::vector<Tensor<2,dim> > > &dduh,
												   const std::vector<Point<dim> >                  &normals,
												   const std::vector<Point<dim> >                  &evaluation_points,
												   const unsigned int							   &material_id,
												   std::vector<Vector<double> >                    &computed_quantities) const;

	/// Returns a vector of strings for the names of the output quantities
	virtual std::vector<std::string> get_names() const;
	/// Method indicating the component interpretation (i.e. is the output part of a vector or is it scalar)
	virtual std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation () const;
	/// Just tells DataOut what to update, in this case the hessians were left off since they are not needed
	virtual UpdateFlags get_needed_update_flags() const { return update_values | update_gradients | update_q_points; };

private:
	/// Copy of materials vector (should be made a pointer)
	std::vector< Material<dim> > mats;
};

template <int dim>
std::vector<std::string> TEPostProcessor<dim>::get_names() const
{
	std::vector<std::string> solution_names;
	switch (dim)
	{
	case 2:
		solution_names.push_back ("T");
		solution_names.push_back ("u_x");
		solution_names.push_back ("u_y");
		solution_names.push_back ("eps_xx");
		solution_names.push_back ("eps_yy");
		solution_names.push_back ("eps_xy");
		solution_names.push_back ("sig_xx");
		solution_names.push_back ("sig_yy");
		solution_names.push_back ("sig_xy");
		break;
	case 3:
		solution_names.push_back ("T");
		solution_names.push_back ("u_x");
		solution_names.push_back ("u_y");
		solution_names.push_back ("u_z");
		solution_names.push_back ("eps_xx");
		solution_names.push_back ("eps_yy");
		solution_names.push_back ("eps_zz");
		solution_names.push_back ("eps_yz");
		solution_names.push_back ("eps_xz");
		solution_names.push_back ("eps_xy");
		solution_names.push_back ("sig_xx");
		solution_names.push_back ("sig_yy");
		solution_names.push_back ("sig_zz");
		solution_names.push_back ("sig_yz");
		solution_names.push_back ("sig_xz");
		solution_names.push_back ("sig_xy");
		break;
	default:
		Assert (false, ExcNotImplemented());
		break;
	}

	return solution_names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> TEPostProcessor<dim>::get_data_component_interpretation () const {
	std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
	// For now set all of the components to scalar, but change later?
	if (dim == 2) {
		for (unsigned int i = 0; i < 1+dim+6; i++)
			interpretation.push_back (DataComponentInterpretation::component_is_scalar);
	}
	else if (dim ==3) {
		for (unsigned int i = 0; i < 1+dim+12; i++)
			interpretation.push_back (DataComponentInterpretation::component_is_scalar);
	}
	return interpretation;
}

template <int dim>
void TEPostProcessor<dim>::compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
															 const std::vector<std::vector<Tensor<1,dim> > > &duh,
															 const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
															 const std::vector<Point<dim> >                  &/*normals*/,
															 const std::vector<Point<dim> >                  &/*evaluation_points*/,
															 std::vector<Vector<double> >                    &computed_quantities) const
{
	// This should not be used now
	Assert (false, ExcMessage("The TEPostProcessor was used with the DataOut class, be sure to use the TEDataOut class"));
}

template <int dim>
void TEPostProcessor<dim>::compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
															 const std::vector<std::vector<Tensor<1,dim> > > &duh,
															 const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
															 const std::vector<Point<dim> >                  &/*normals*/,
															 const std::vector<Point<dim> >                  &/*evaluation_points*/,
															 const unsigned int								 &material_id,
															 std::vector<Vector<double> >                    &computed_quantities) const
{
	const unsigned int n_quadrature_points = uh.size();
	Assert (duh.size() == n_quadrature_points, ExcInternalError());
	Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
	Assert (uh[0].size() == dim+1, ExcInternalError());

	for (unsigned int q=0; q<n_quadrature_points; ++q) {
		// Variable used to keep track of the last used index
		unsigned int curr_ind = 0;

		// Save temperature (index 0)
		double temperature = uh[q](0);
		computed_quantities[q](curr_ind) = temperature;
		curr_ind += 1;

		// Save displacement and compute grad_u
		// Displacement is indices 1-dim
		Tensor<2, dim> grad_u;
		for (unsigned int i = 0; i < dim; i++) {
			computed_quantities[q](i+curr_ind) = uh[q](i+curr_ind);
			grad_u[i] = duh[q][i+curr_ind];
		}
		curr_ind += dim;

		// Compute strain, symmetrize returns (grad_u + transpose(grad_u))/2
		const SymmetricTensor<2,dim> total_strain = symmetrize (grad_u);
		// Output in Voigt notation
		// For 3D its: (indices dim+1 - dim+6)
		// eps_xx eps_yy eps_zz eps_yz eps_xz eps_xy
		// For 2D its: (indices dim+1 - dim+3)
		// eps_xx eps_yy eps_xy
		// NOTE: There has to be a more systematic way of doing this for the future, reevaluate later
		if (dim == 2) {
			computed_quantities[q](curr_ind) = total_strain[0][0];
			computed_quantities[q](curr_ind+1) = total_strain[1][1];
			computed_quantities[q](curr_ind+2) = total_strain[0][1];
			curr_ind += 3;
		}
		else if (dim == 3) {
			computed_quantities[q](curr_ind) = total_strain[0][0];
			computed_quantities[q](curr_ind+1) = total_strain[1][1];
			computed_quantities[q](curr_ind+2) = total_strain[2][2];
			computed_quantities[q](curr_ind+3) = total_strain[1][2];
			computed_quantities[q](curr_ind+4) = total_strain[0][2];
			computed_quantities[q](curr_ind+5) = total_strain[0][1];
			curr_ind += 6;
		}

		// Get cell's material properties
		SymmetricTensor<4, dim> stiffness;
		SymmetricTensor<2, dim> alpha;

		bool mat_found = false;
		for (unsigned int mat_ind = 0; mat_ind < mats.size(); mat_ind++) {
			if (mats[mat_ind].get_id() == material_id) {
				stiffness = mats[mat_ind].get_stiffness();
				alpha = mats[mat_ind].get_expansion();
				mat_found = true;
				break;
			}
		}
		Assert(mat_found, ExcMessage("Material not found in TEPostProcessor."))

		// Compute eigenstrain
		// How do I get the coefficient of thermal expansion matrix?
		// And get the stiffness matrix?
		const SymmetricTensor<2,dim> elastic_strain = total_strain-alpha*temperature;
		// Compute stress
		const SymmetricTensor<2,dim> stress = stiffness * elastic_strain;
		if (dim == 2) {
			computed_quantities[q](curr_ind) = stress[0][0];
			computed_quantities[q](curr_ind+1) = stress[1][1];
			computed_quantities[q](curr_ind+2) = stress[0][1];
			curr_ind += 3;
		}
		else if (dim == 3) {
			computed_quantities[q](curr_ind) = stress[0][0];
			computed_quantities[q](curr_ind+1) = stress[1][1];
			computed_quantities[q](curr_ind+2) = stress[2][2];
			computed_quantities[q](curr_ind+3) = stress[1][2];
			computed_quantities[q](curr_ind+4) = stress[0][2];
			computed_quantities[q](curr_ind+5) = stress[0][1];
			curr_ind += 6;
		}

	}
}

}

#endif /* TEPOSTPROCESSOR_H_ */
