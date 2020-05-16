#include "../enum_and_typedefs.h"
#include <stdio.h>

int curl_dual(Dual_vector_field in_field, Vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int layer_index, h_index, index_0, index_1, index_2, index_3, index_4, index_5, sign_0, sign_1, sign_2, sign_3, sign_4, sign_5;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if(h_index < NUMBER_OF_VECTORS_V)
        {
            index_0 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> vorticity_indices[6*h_index + 0];
            index_1 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> vorticity_indices[6*h_index + 1];
            index_2 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> vorticity_indices[6*h_index + 2];
            index_3 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> vorticity_indices[6*h_index + 3];
            index_4 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> vorticity_indices[6*h_index + 4];
            index_5 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> vorticity_indices[6*h_index + 5];
            sign_0 = grid -> vorticity_signs[6*h_index + 0];
            sign_1 = grid -> vorticity_signs[6*h_index + 1];
            sign_2 = grid -> vorticity_signs[6*h_index + 2];
            sign_3 = grid -> vorticity_signs[6*h_index + 3];
            sign_4 = grid -> vorticity_signs[6*h_index + 4];
            sign_5 = grid -> vorticity_signs[6*h_index + 5];
            out_field[i] = 1/grid -> area[i]*(dualgrid -> normal_distance[index_0]*sign_0*in_field[index_0] + dualgrid -> normal_distance[index_1]*sign_1*in_field[index_1] + dualgrid -> normal_distance[index_2]*sign_2*in_field[index_2] + dualgrid -> normal_distance[index_3]*sign_3*in_field[index_3] + dualgrid -> normal_distance[index_4]*sign_4*in_field[index_4] + dualgrid -> normal_distance[index_5]*sign_5*in_field[index_5]);
        }
        else
        {
            index_0 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> h_curl_indices[4*(h_index - NUMBER_OF_VECTORS_V) + 0];
            index_1 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> h_curl_indices[4*(h_index - NUMBER_OF_VECTORS_V) + 1];
            index_2 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> h_curl_indices[4*(h_index - NUMBER_OF_VECTORS_V) + 2];
            index_3 = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + grid -> h_curl_indices[4*(h_index - NUMBER_OF_VECTORS_V) + 3];
            sign_0 = grid -> h_curl_signs[4*(h_index - NUMBER_OF_VECTORS_V) + 0];
            sign_1 = grid -> h_curl_signs[4*(h_index - NUMBER_OF_VECTORS_V) + 1];
            sign_2 = grid -> h_curl_signs[4*(h_index - NUMBER_OF_VECTORS_V) + 2];
            sign_3 = grid -> h_curl_signs[4*(h_index - NUMBER_OF_VECTORS_V) + 3];
            out_field[i] = 1/grid -> area[i]*(dualgrid -> normal_distance[index_0]*sign_0*in_field[index_0] + dualgrid -> normal_distance[index_1]*sign_1*in_field[index_1] + dualgrid -> normal_distance[index_2]*sign_2*in_field[index_2] + dualgrid -> normal_distance[index_3]*sign_3*in_field[index_3]);
        }
    }
    return 0;
}
