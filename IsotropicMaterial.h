/*
 * IsotropicMaterial.h
 *
 *  Created on: Mar 14, 2013
 *      Author: kballard
 */

#ifndef ISOTROPICMATERIAL_H_
#define ISOTROPICMATERIAL_H_

namespace FEASolverNS
{

class IsotropicMaterial {
public:
	IsotropicMaterial(unsigned int material_id, double lambda, double mu) : mat_id(material_id), lambda_val(lambda), mu_val(mu) {};

	unsigned int get_id() { return mat_id; }
	double get_lambda() { return lambda_val; }
	double get_mu() { return mu_val; }

private:
	unsigned int mat_id;
	double lambda_val, mu_val;
};

}

#endif /* ISOTROPICMATERIAL_H_ */
