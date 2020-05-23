#include "../enum_and_typedefs.h"
#include <stdio.h>

int grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double vertical_gradient;
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
    for (int i = (NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)*NUMBER_OF_VECTORS_H; i < NUMBER_OF_H_VECTORS; ++i)
    {
    	layer_index = i/NUMBER_OF_VECTORS_H;
    	h_index = i - layer_index*NUMBER_OF_VECTORS_H;
    	vertical_gradient = 0.25*(out_field[grid -> to_index[h_index] + layer_index*NUMBER_OF_VECTORS_PER_LAYER] + out_field[grid -> to_index[h_index] + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER] + out_field[grid -> from_index[h_index] + layer_index*NUMBER_OF_VECTORS_PER_LAYER] + out_field[grid -> from_index[h_index] + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]);
    	out_field[i] = out_field[i] - (grid -> z_scalar[grid -> to_index[h_index] + layer_index*NUMBER_OF_SCALARS_H] - grid -> z_scalar[grid -> from_index[h_index] + layer_index*NUMBER_OF_SCALARS_H])/grid -> normal_distance[i]*vertical_gradient;
    	
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
        out_field[i] = out_field[i + NUMBER_OF_VECTORS_PER_LAYER];
    for (int i = NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V; i < NUMBER_OF_VECTORS; ++i)
        out_field[i] = out_field[i - NUMBER_OF_VECTORS_PER_LAYER];
    return 0;
}
