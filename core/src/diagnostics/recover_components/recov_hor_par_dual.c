#include "../../enum_and_typedefs.h"
#include <stdio.h>

int recov_hor_par_dual(Dual_vector_field in_field, long layer_index, long h_index, double *component, Grid *grid)
{
    *component = 0;
    for (int i = 0; i < 2; ++i)
        *component += grid -> recov_hor_par_dual_weight[2*h_index + i]*in_field[layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> recov_hor_par_dual_index[2*h_index + i]];
    return 0;
}