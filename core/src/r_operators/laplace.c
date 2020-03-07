#include "../enum_and_typedefs.h"
#include "r_operators.h"

void laplace(Scalar_field in_field, Scalar_field out_field)
{
    Vector_field between;
    grad(in_field, between);
    divergence(between, out_field);
}
