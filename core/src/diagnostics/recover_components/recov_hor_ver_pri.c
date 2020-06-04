#include "../../enum_and_typedefs.h"

int recov_hor_ver_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    for (int i = 0; i < 4; ++i)
        *component += grid -> recov_hor_ver_pri_weight[4*h_index + i]*in_field[layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + i]];
    return 0;
}
