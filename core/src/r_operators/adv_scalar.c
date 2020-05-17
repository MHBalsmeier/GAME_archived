#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>

int adv_scalar(Scalar_field property, Vector_field velocity, Scalar_field result, Grid *grid, Dualgrid *dualgrid)
{
    Vector_field *gradient = malloc(sizeof(Vector_field));
    int retval = grad(property, *gradient, grid);
    retval = inner(velocity, *gradient, result, grid, dualgrid);
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
        result[i] = -result[i];
    free(gradient);
    return retval;
}
