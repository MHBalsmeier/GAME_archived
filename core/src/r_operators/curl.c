#include "../enum_and_typedefs.h"
#include <stdio.h>

int curl(Vector_field in_field, Dual_vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    long layer_index, h_index;
    long index_0, index_1, index_2, index_3;
    short sign_0, sign_1, sign_2, sign_3;
    double distance_0, distance_1, distance_2, distance_3, rhombus_circ, rhombus_area, check;
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
        {
            out_field[i] = 0;
            for (int j = 0; j < 3; ++j)
            {
                index_0 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + 0];
                index_1 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + 1];
                index_2 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + 2];
                sign_0 = dualgrid -> vorticity_signs[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + 0];
                sign_1 = dualgrid -> vorticity_signs[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + 1];
                sign_2 = dualgrid -> vorticity_signs[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + 2];
                rhombus_circ = grid -> normal_distance[index_0]*sign_0*in_field[index_0] + grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_2]*sign_2*in_field[index_2];
                check = rhombus_circ;
                index_0 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[3*dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + 0];
                index_1 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[3*dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + 1];
                index_2 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[3*dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + 2];
                sign_0 = dualgrid -> vorticity_signs[3*dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + 0];
                sign_1 = dualgrid -> vorticity_signs[3*dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + 1];
                sign_2 = dualgrid -> vorticity_signs[3*dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + 2];
                rhombus_circ += grid -> normal_distance[index_0]*sign_0*in_field[index_0] + grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_2]*sign_2*in_field[index_2];
                rhombus_area = dualgrid -> area[i] + dualgrid -> area[dualgrid -> adjacent_scalar_indices_dual_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j] + NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
                out_field[i] += 1.0/3.0*rhombus_circ/rhombus_area;
            }
        }
        else
        {
            if (layer_index == 0 || layer_index == NUMBER_OF_LAYERS)
            {
                index_1 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_3 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                out_field[i] = 1/dualgrid -> area[i]*(grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_3]*sign_3*in_field[index_3]);
            }
            else
            {
                index_0 = layer_index*NUMBER_OF_VECTORS_PER_LAYER - NUMBER_OF_VECTORS_H + dualgrid -> h_curl_indices[4*h_index + 0];
                index_1 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_2 = layer_index*NUMBER_OF_VECTORS_PER_LAYER - NUMBER_OF_VECTORS_H + dualgrid -> h_curl_indices[4*h_index + 2];
                index_3 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_0 = dualgrid -> h_curl_signs[4*h_index + 0];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_2 = dualgrid -> h_curl_signs[4*h_index + 2];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                out_field[i] = 1/dualgrid -> area[i]*(grid -> normal_distance[index_0]*sign_0*in_field[index_0] + grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_2]*sign_2*in_field[index_2] + grid -> normal_distance[index_3]*sign_3*in_field[index_3]);
            }
        }
    }
    return 0;
}
