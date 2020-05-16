#include "../enum_and_typedefs.h"
#include <stdio.h>

int divergence(Vector_field in_field, Scalar_field out_field, Grid *grid, int allow_surface_flux)
{
    int number_of_edges, layer_index, h_index;
    double comp_h, comp_v;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        comp_h = 0;
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        number_of_edges = 6;
        if (h_index < NUMBER_OF_PENTAGONS)
            number_of_edges = 5;
        for (int j = 0; j < number_of_edges; ++j)
            comp_h += in_field[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]*grid -> adjacent_signs_h[6*h_index + j]*grid -> area[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
        if (layer_index == 0)
            comp_v = -in_field[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]*grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        else if (layer_index == NUMBER_OF_LAYERS - 1)
            comp_v = in_field[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER]*grid -> area[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER] - allow_surface_flux*in_field[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]*grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        else
            comp_v = in_field[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER]*grid -> area[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER] - in_field[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]*grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        out_field[i] = 1/grid -> volume[i]*(comp_h + comp_v);
    }
    return 0;
}
