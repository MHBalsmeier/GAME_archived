#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>

int laplace_vec(Vector_field in_field, Vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    Scalar_field *between_0 = malloc(sizeof(Scalar_field));
    divergence(in_field, *between_0, grid);
    Vector_field *between_1 = malloc(sizeof(Vector_field));
    grad(*between_0, *between_1, grid);
    free(between_0);
    Dual_vector_field *between_2 = malloc(sizeof(Dual_vector_field));
    rot(in_field, *between_2, grid, dualgrid);
    Vector_field *between_3 = malloc(sizeof(Vector_field));
    rot_dual(*between_2, *between_3, grid, dualgrid);
    free(between_2);
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        out_field[i] = (*between_1)[i] - (*between_3)[i];
    free(between_1);
    free(between_3);
    return 0;
}
