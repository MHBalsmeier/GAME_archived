#include "../enum_and_typedefs.h"
#include "r_operators.h"

void adv_s(Vector_field wind, Scalar_field prop, Scalar_field out_field)
{
	Scalar_field advection;
	Vector_field gradient_field;
	grad(prop, gradient_field);
	scalar_product(wind, gradient_field, advection);
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
		out_field[i] = -advection[i];
}
