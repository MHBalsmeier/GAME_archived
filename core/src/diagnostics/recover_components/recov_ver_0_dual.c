#include "../../enum_and_typedefs.h"

int recov_ver_0_dual(Dual_vector_field in_field, long layer_index, long h_index, double *component, Grid *grid)
{
    *component = 0;
    for (int i = 0; i < 6; ++i)
        *component += grid -> recov_ver_0_dual_weight[6*h_index + i]*in_field[layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> recov_ver_0_dual_index[6*h_index + i]];
    return 0;
}
