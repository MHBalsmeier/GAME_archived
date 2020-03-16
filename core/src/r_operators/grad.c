#include "../enum_and_typedefs.h"
#include <stdio.h>

void grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
    long layer_index, h_index, lower_index, upper_index;
    for (int i = NUMBER_OF_VECTORS_V; i < NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
            out_field[i] = (in_field[grid -> to_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H] - in_field[grid -> from_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H])/grid -> normal_distance[i];
        else
        {
            lower_index = h_index + layer_index*NUMBER_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H;
            out_field[i] = (in_field[upper_index] - in_field[lower_index])/grid -> normal_distance[i];
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
        out_field[i] = out_field[i + NUMBER_OF_VECTORS_PER_LAYER];
    for (int i = NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V; i < NUMBER_OF_VECTORS; ++i)
        out_field[i] = out_field[i - NUMBER_OF_VECTORS_PER_LAYER];
}
