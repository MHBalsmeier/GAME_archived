/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include <stdio.h>

int trsk_modified(Vector_field in_field_0, Curl_field in_field_1, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    for (int i = 0; i < 10; ++i)
        *component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NUMBER_OF_DUAL_VECTORS_H + layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]];
    return 0;
}
