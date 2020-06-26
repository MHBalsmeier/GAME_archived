/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"

int recov_ver_0_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    if (layer_index > 0 && layer_index < NUMBER_OF_LAYERS)
    {
        for (int i = 0; i < 6; ++i)
        {
            *component += 0.5*grid -> recov_ver_0_pri_weight[6*h_index + i]*in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
            *component += 0.5*grid -> recov_ver_0_pri_weight[6*h_index + i]*in_field[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
        }
    }
    if (layer_index == 0)
    {
        for (int i = 0; i < 6; ++i)
            *component += grid -> recov_ver_0_pri_weight[6*h_index + i]*in_field[NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
    }
    if (layer_index == NUMBER_OF_LAYERS)
    {
        for (int i = 0; i < 6; ++i)
            *component += grid -> recov_ver_0_pri_weight[6*h_index + i]*in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
    }
    return 0;
}
