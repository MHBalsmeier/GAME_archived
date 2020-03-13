#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>
#include <stdio.h>

int adv_s(Vector_field wind, Scalar_field prop, Scalar_field out_field, Grid *grid)
{
    Scalar_field *minus_advection = malloc(sizeof(Scalar_field));
    Vector_field *gradient_field = malloc(sizeof(Vector_field));
    grad(prop, *gradient_field, grid);
    scalar_product(wind, *gradient_field, *minus_advection, grid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        out_field[i] = -(*minus_advection)[i];
    free(minus_advection);
    return 0;
}
