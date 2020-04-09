#include "../enum_and_typedefs.h"
#include <stdio.h>

int inner(Vector_field in_field_0, Vector_field in_field_1, Scalar_field out_field, Grid *grid)
{
    long layer_index, h_index;
    double comp_h, comp_v, x_0, x_1, y_0, y_1;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        x_0 = 0;
        x_1 = 0;
        y_0 = 0;
        y_1 = 0;
        for (int i = 0; i < 6; ++i)
        {
            x_0 += grid -> recov_ver_0_pri_weight[6*h_index + i]*in_field_0[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_0_pri_index[6*h_index + i]];
            x_1 += grid -> recov_ver_0_pri_weight[6*h_index + i]*in_field_1[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_0_pri_index[6*h_index + i]];
            y_0 += grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field_0[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_1_pri_index[6*h_index + i]];
            y_1 += grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field_1[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_1_pri_index[6*h_index + i]];
        }
        comp_h = x_0*x_1 + y_0*y_1;
        comp_v = in_field_0[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER]*in_field_1[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
        comp_v += in_field_0[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]*in_field_1[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        out_field[i] = comp_h + 0.5*comp_v;
    }
    return 0;
}
