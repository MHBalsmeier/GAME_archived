/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"

int recov_ver_1_curl(Curl_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    for (int i = 0; i < 6; ++i)
        *component += grid -> recov_ver_1_curl_weight[6*h_index + i]*in_field[layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]];
    return 0;
}
