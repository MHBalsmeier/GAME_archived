#include "../enum_and_typedefs.h"
#include <stdio.h>

void divergence(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
    long layer_index, h_index;
    double comp_h, comp_v;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        if (h_index < NUMBER_OF_PENTAGONS)
        {
            for (int j = 0; j < 5; ++j)
                comp_h += in_field[layer_index*NUMBER_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]]*grid -> adjacent_signs_h[6*h_index + j]*grid -> area[layer_index*NUMBER_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
        }
        else
        {
            for (int j = 0; j < 6; ++j)
                comp_h += in_field[layer_index*NUMBER_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]]*grid -> adjacent_signs_h[6*h_index + j]*grid -> area[layer_index*NUMBER_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
        }
        comp_v = in_field[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H]*grid -> area[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H] - in_field[h_index + layer_index*NUMBER_OF_SCALARS_H]*grid -> area[h_index + layer_index*NUMBER_OF_SCALARS_H];
        out_field[i] = (1/grid -> volume[i])*(comp_h + comp_v);
    }
}
