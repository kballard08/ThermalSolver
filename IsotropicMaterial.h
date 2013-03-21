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
	IsotropicMaterial(unsigned int material_id, double lambda, double mu, double alpha);

	unsigned int get_id() { return mat_id; }
	double get_lambda() { return lambda_val; }
	double get_mu() { return mu_val; }
	Tensor<2, 3> get_alpha() { return alpha_ten; }

private:
	unsigned int mat_id;
	double lambda_val, mu_val;
	// Rank 2 tensor 3 indicies (i.e. a 3x3 matrix)
	Tensor<2, 3> alpha_ten;
};

IsotropicMaterial::IsotropicMaterial(unsigned int material_id, double lambda, double mu, double alpha) : mat_id(material_id), lambda_val(lambda), mu_val(mu), alpha_ten()
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i ==j)
				alpha_ten[i][j] = alpha;
			else
				alpha_ten[i][j] = 0.0;
		}
	}
}

}

#endif /* ISOTROPICMATERIAL_H_ */