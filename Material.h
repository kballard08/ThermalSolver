/*
 * Material.h
 *
 *  Created on: Apr 22, 2013
 *      Author: kballard
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

namespace FEASolverNS
{

template<int dim>
class Material
{
public:
	// Virtual constructor and destructor
	Material(const unsigned int &id) : mat_id(id) {};
	virtual ~Material() {};

	// Virtual methods that will be called by anything that uses a material
	// Gets material id
	virtual unsigned int get_id() { return mat_id; };
	// Gets the 4th order stiffness tensor
	virtual SymmetricTensor<4, dim> get_stiffness() { return stiffness_tensor; };
	// Gets the 2nd order coefficient of thermal expansion tensor
	virtual SymmetricTensor<2, dim> get_expansion() { return alpha_tensor; };
	// Gets the 2nd order coefficient of thermal conduction tensor
	virtual Tensor<2, dim> get_conduction() { return k_tensor; };

protected:
	// The material id that will correspond to the material id in the mesh
	unsigned int mat_id;
	// Rank 2 tensor of size dim (conduction)
	// NOTE: This tensor is not stored as a symmetric for ease of use during assembly
	Tensor<2, dim> k_tensor;
	// Rank 2 tensor of size dim (expansion)
	SymmetricTensor<2, dim> alpha_tensor;
	// Rank 4 tensor of size dim (stiffness)
	SymmetricTensor<4, dim> stiffness_tensor;
}; // Material class

} // FEASolverNS

#endif /* MATERIAL_H_ */
