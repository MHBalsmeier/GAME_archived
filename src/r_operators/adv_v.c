#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"

void adv_v(Vector_field wind, Vector_field out_field)
{
	Scalar_field kinetic_energy;
	scalar_product(wind,wind,kinetic_energy);
	Vector_field grad_kinetic_energy;
	grad(kinetic_energy,grad_kinetic_energy);
	Dual_vector_field rot_wind;
	rot(wind,rot_wind);
	Vector_field product_curl;
	vector_product(wind,rot_wind,product_curl);
	for (int i = 0; i<=NUMBER_OF_VECTORS-1; ++i)
	{
		out_field[i] = -0.5*grad_kinetic_energy[i] + product_curl[i];
	}
}
