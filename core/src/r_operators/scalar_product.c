#include "../enum_and_typedefs.h"
#include <stdio.h>

int scalar_product(Vector_field in_field_1, Vector_field in_field_2, Scalar_field out_field, Grid *grid)
{
    double comp_h, comp_v, factor_h, factor_v;
    long layer_index, h_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        comp_h, comp_v = 0;
        factor_h = 2.0/6.0;
        factor_v = 0.5;
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        for (int j = 0; j < 6; ++j)
        {
            if (j < 5 || h_index >= NUMBER_OF_PENTAGONS)
            comp_h += in_field_1[layer_index*NUMBER_OF_SCALARS_H + grid -> adjacent_vector_indices_h[6*h_index + j]]*in_field_2[layer_index*NUMBER_OF_SCALARS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
        }
        if (layer_index >= 1)
            comp_v += in_field_1[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H]*in_field_2[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H];
        else
            factor_v = 1;
        if (layer_index < NUMBER_OF_LAYERS - 1)
            comp_v += in_field_1[h_index + layer_index*NUMBER_OF_SCALARS_H]*in_field_2[h_index + layer_index*NUMBER_OF_SCALARS_H];
        else
            factor_v = 1;
        if (h_index < NUMBER_OF_PENTAGONS)
            factor_h = 2.0/5.0;
        out_field[i] = factor_h*comp_h + factor_v*comp_v;
    }   
}
