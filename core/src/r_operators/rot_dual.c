#include "../enum_and_typedefs.h"

int rot_dual(Dual_vector_field in_field, Vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int vert_index, floor_index, h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        vert_index = i/(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        floor_index = vert_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        h_index = i - floor_index;
        if(h_index < NUMBER_OF_SCALARS_H)
        {
            int vertical_index;
            vertical_index = dualgrid -> to_index[i];
            int index_1, index_2, index_3, index_4, index_5, index_6, sign_1, sign_2, sign_3, sign_4, sign_5, sign_6;
            index_1 = grid -> vorticity_indices[6*vertical_index + 0];
            index_2 = grid -> vorticity_indices[6*vertical_index + 1];
            index_3 = grid -> vorticity_indices[6*vertical_index + 2];
            index_4 = grid -> vorticity_indices[6*vertical_index + 3];
            index_5 = grid -> vorticity_indices[6*vertical_index + 4];
            sign_1 = grid -> vorticity_signs[6*vertical_index + 0];
            sign_2 = grid -> vorticity_signs[6*vertical_index + 1];
            sign_3 = grid -> vorticity_signs[6*vertical_index + 2];
            sign_4 = grid -> vorticity_signs[6*vertical_index + 3];
            sign_5 = grid -> vorticity_signs[6*vertical_index + 4];
            if(h_index > NUMBER_OF_PENTAGONS)
            {
                out_field[i] = (1/grid -> area[i])*(dualgrid -> normal_distance[index_1]*sign_1*in_field[index_1] + dualgrid -> normal_distance[index_2]*sign_2*in_field[index_2] + dualgrid -> normal_distance[index_3]*sign_3*in_field[index_3] + dualgrid -> normal_distance[index_4]*sign_4*in_field[index_4] + dualgrid -> normal_distance[index_5]*sign_5*in_field[index_5]);
            }
            else
            {
                index_6 = grid -> vorticity_indices[6*vertical_index + 5];
                sign_6 = grid -> vorticity_signs[6*vertical_index + 5];
                out_field[i] = (1/grid -> area[i])*(dualgrid -> normal_distance[index_1]*sign_1*in_field[index_1] + dualgrid -> normal_distance[index_2]*sign_2*in_field[index_2] + dualgrid -> normal_distance[index_3]*sign_3*in_field[index_3] + dualgrid -> normal_distance[index_4]*sign_4*in_field[index_4] + dualgrid -> normal_distance[index_5]*sign_5*in_field[index_5] + dualgrid -> normal_distance[index_6]*sign_6*in_field[index_6]);
            }
        }
        else
        {
            int horizontal_index;
            horizontal_index = grid -> to_index[i];
            int index_1, index_2, index_3, index_4, sign_1, sign_2, sign_3, sign_4;
            double distance_1, distance_2, distance_3, distance_4;
            index_1 = grid -> h_curl_indices[4*horizontal_index + 0];
            index_2 = grid -> h_curl_indices[4*horizontal_index + 1];
            index_3 = grid -> h_curl_indices[4*horizontal_index + 2];
            index_4 = grid -> h_curl_indices[4*horizontal_index + 3];
            sign_1 = dualgrid -> h_curl_signs[4*horizontal_index + 0];
            sign_2 = dualgrid -> h_curl_signs[4*horizontal_index + 1];
            sign_3 = dualgrid -> h_curl_signs[4*horizontal_index + 2];
            sign_4 = dualgrid -> h_curl_signs[4*horizontal_index + 3];
            distance_1 = dualgrid -> normal_distance[index_1];
            distance_2 = dualgrid -> normal_distance[index_2];
            distance_3 = dualgrid -> normal_distance[index_3];
            distance_4 = dualgrid -> normal_distance[index_4];
            out_field[i] = (1/grid -> area[i])*(distance_1*sign_1*in_field[index_1] + distance_2*sign_2*in_field[index_2] + distance_3*sign_3*in_field[index_3] + distance_4*sign_4*in_field[index_4]);
        }
    }
    return 0;
}
