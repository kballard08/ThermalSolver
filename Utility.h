/*
 * Utility.h
 *
 *  Created on: Mar 1, 2013
 *      Author: kballard
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>

#include <iostream>
#include <string>

namespace FEASolverNS
{

enum Verbosity {
	MIN_V,
	MAX_V,
	DEBUG_V
};

void Status(std::string message, int verbosity, int verbosity_required) {
	if (verbosity >= verbosity_required)
		std::cout << message << std::endl;
}

double kron_delta(const unsigned int &i, const unsigned int &j) {
	if (i == j)
		return 1.0;
	return 0.0;
}

double kron_delta(const int &i, const int &j) {
	if (i == j)
		return 1.0;
	return 0.0;
}

// * operator for Vector and vector
double operator * (const dealii::Vector<double> &v1, const std::vector<double> &v2) {
	// TODO: Add assertion here for size
	// Compute dot project
	double sum = 0;
	for (unsigned int i = 0; i < v1.size(); i++)
		sum += v1[i]*v2[i];
	return sum;
}
double operator * (const std::vector<double> &v1, const dealii::Vector<double> &v2) {
	// Just use the other operator for ease of maintainence
	return v2*v1;
}

// * operator for Tensor (rank 1) and Vector
template <int dim>
double operator * (const dealii::Vector<double> &v1, const  dealii::Tensor<1, dim> &v2) {
	// TODO: Add assertion here for size
	// Compute dot project
	double sum = 0;
	for (unsigned int i = 0; i < dim; i++)
		sum += v1[i]*v2[i];
	return sum;
}
template <int dim>
double operator * (const dealii::Tensor<1, dim> &v1, const dealii::Vector<double> &v2) {
	// Just use the other operator for ease of maintainence
	return v2*v1;
}

// * operator for Tensor (rank 4) and Tensor (rank 2)
template <int dim>
dealii::Tensor<2, dim> operator * (const dealii::Tensor<4, dim> &v1, const  dealii::Tensor<2, dim> &v2) {
	// TODO: Add assertion here for size

	dealii::Tensor<2, dim> result;
	// result_ij = v1_ijkl*v2_kl
	for (unsigned int i = 0; i < dim; i++)
		for (unsigned int j = 0; j < dim; j++)
			for (unsigned int k = 0; k < dim; k++)
				for (unsigned int l = 0; l < dim; l++)
					result[i][j] = v1[i][j][k][l]*v2[k][l];
	return result;
}

}

#endif /* UTILITY_H_ */
