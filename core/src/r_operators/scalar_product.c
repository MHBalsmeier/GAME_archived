#include "../enum_and_typedefs.h"

void scalar_product(Vector_field in_field_1, Vector_field in_field_2, Scalar_field out_field)
{
    extern Grid grid;
    double comp_h, comp_v, factor_h, factor_v;
    long layer_index, horizontal_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        comp_h, comp_v = 0;
        factor_h = 2.0/6.0;
        factor_v = 0.5;
        layer_index = i/NUMBER_OF_SCALARS_H;
        horizontal_index = i - layer_index*NUMBER_OF_SCALARS_H;
        for (int j = 0; j < 6; j++)
        {
            if (j < 5 || horizontal_index >= NUMBER_OF_PENTAGONS)
            comp_h += in_field_1[grid.adjacent_vector_indices_h[6*i + j]]*in_field_2[grid.adjacent_vector_indices_h[6*i + j]];
        }
        if (layer_index >= 1)
            comp_v += in_field_1[grid.adjacent_vector_index_upper[i]]*in_field_2[grid.adjacent_vector_index_upper[i]];
        else
            factor_v = 1;
        if (layer_index < NUMBER_OF_LAYERS - 1)
            comp_v += in_field_1[grid.adjacent_vector_index_lower[i]]*in_field_2[grid.adjacent_vector_index_lower[i]];
        else
            factor_v = 1;
        if (horizontal_index < NUMBER_OF_PENTAGONS)
            factor_h = 2.0/5.0;
        out_field[i] = factor_h*comp_h + factor_v*comp_v;
    }   
}
