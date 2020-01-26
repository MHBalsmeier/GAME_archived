#include "../enum_and_typedefs.h"
#include "r_operators.h"

void laplace_vec(Vector_field in_field, Vector_field out_field)
{
	Scalar_field between;
	divergence(in_field, between);
	Vector_field between_1;
	grad(between, between_1);
	Dual_vector_field between_2;
	rot(in_field, between_2);
	Vector_field between_3;
	rot_dual(between_2, between_3);
	for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
		out_field[i] = between_1[i] - between_3[i];
}
