/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int scalar_times_vector(Scalar_field in_field_0, Vector_field in_field_1, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double scalar_value;
    for (int i = NUMBER_OF_VECTORS_V; i < NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
            scalar_value = 0.5*(in_field_0[grid -> to_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H] + in_field_0[grid -> from_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H]);
        else
        {
            lower_index = h_index + layer_index*NUMBER_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H;
            scalar_value = 0.5*(in_field_0[upper_index] + in_field_0[lower_index]);
        }
        out_field[i] = scalar_value*in_field_1[i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        scalar_value = in_field_0[i];
        out_field[i] = scalar_value*in_field_1[i];
    }
    for (int i = NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H;
        scalar_value = in_field_0[upper_index];
        out_field[i] = scalar_value*in_field_1[i];
    }
    return 0;
}
