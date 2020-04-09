#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>

int laplace(Scalar_field in_field, Scalar_field out_field, Grid *grid)
{
    Vector_field *between = malloc(sizeof(Vector_field));
    grad(in_field, *between, grid);
    divergence(*between, out_field, grid, 0);
    free(between);
    return 0;
}
