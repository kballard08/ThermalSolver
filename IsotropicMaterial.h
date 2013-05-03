/*
 * IsotropicMaterial.h
 *
 *  Created on: Mar 14, 2013
 *      Author: kballard
 */

#ifndef ISOTROPICMATERIAL_H_
#define ISOTROPICMATERIAL_H_

#include "Material.h"
#include "Utility.h"

namespace FEASolverNS
{

/**
 * Class derived from Material class describing an isotropic material.  It is very easy to add transversely isotropic and
 * orthotropic material classes as needed.
 */
template<int dim>
class IsotropicMaterial : public Material<dim>
{
public:
	/// Constructor
	IsotropicMaterial(const unsigned int &material_id, const double &E, const double &nu, const double &alpha, const double &k);
	/// Destructor
	virtual ~IsotropicMaterial() {};

private:
	/// Isotropic constants (I kept these around for debug purposes, but these can be moved to just in the constructor now)
	double lambda, mu;
};

template<int dim>
IsotropicMaterial<dim>::IsotropicMaterial(const unsigned int &material_id, const double &E, const double &nu, const double &alpha, const double &k)
	: Material<dim>(material_id)
{
	// Lambda and mu material constants calculated form the relations assuming Hooke's law
	// Obviously only valid for isotropic materials.
	lambda = E*nu/((1+nu)*(1-2*nu));
	mu = E/(2*(1+nu));

	// Form coefficient of thermal conduction tensor
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			Material<dim>::k_tensor[i][j] = k*kron_delta(i,j);

	// Form coefficient of thermal expansion tensor
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			Material<dim>::alpha_tensor[i][j] = alpha*kron_delta(i,j);

	// Form stiffness tensor
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			for (int k = 0; k < dim; k++)
				for (int l = 0; l < dim; l++)
					Material<dim>::stiffness_tensor[i][j][k][l] = lambda*kron_delta(i,j)*kron_delta(k,l) + mu*(kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k));
}

}

#endif /* ISOTROPICMATERIAL_H_ */
