#include "../../enum_and_typedefs.h"

int recov_ver_1_pri(Vector_field in_field, long layer_index, long h_index, double *component, Grid *grid)
{
    *component = 0;
    if (layer_index > 0 && layer_index < NUMBER_OF_LAYERS)
    {
        for (int i = 0; i < 12; ++i)
            *component += grid -> recov_ver_1_pri_weight[12*h_index + i]*in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_1_pri_index[12*h_index + i]];
    }
    if (layer_index == 0)
    {
        for (int i = 0; i < 6; ++i)
            *component += 2*grid -> recov_ver_1_pri_weight[12*h_index + i]*in_field[NUMBER_OF_VECTORS_V + grid -> recov_ver_1_pri_index[12*h_index + i]];
    }
    if (layer_index == NUMBER_OF_LAYERS)
    {
        for (int i = 0; i < 6; ++i)
            *component += 2*grid -> recov_ver_1_pri_weight[12*h_index + i]*in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_1_pri_index[12*h_index + i]];
    }
    return 0;
}

