#include "../../enum_and_typedefs.h"

int recov_hor_par_pri(Vector_field in_field, long layer_index, long h_index, double *component, Grid *grid)
{
    *component = 0;
    for (int i = 0; i < 4; ++i)
        *component += grid -> recov_hor_par_pri_weight[4*h_index + i]*in_field[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> recov_hor_par_pri_index[4*h_index + i]];
    return 0;
}
