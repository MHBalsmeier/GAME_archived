#include "../enum_and_typedefs.h"

void rot(Vector_field in_field, Dual_vector_field out_field)
{
    extern Grid grid;
	extern Dualgrid dualgrid;
    int vertical_index, layer_index, h_index;
	for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
	{
        layer_index = i/(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_SCALARS_H);
        h_index = i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_SCALARS);
        vertical_index = 0;
		if(h_index < NUMBER_OF_DUAL_VECTORS_H)
		{
			int index_1, index_2, index_3, sign_1, sign_2, sign_3;
			index_1 = dualgrid.vorticity_indices[3*vertical_index + 0];
			index_2 = dualgrid.vorticity_indices[3*vertical_index + 1];
			index_3 = dualgrid.vorticity_indices[3*vertical_index + 2];
			sign_1 = dualgrid.vorticity_signs[3*vertical_index + 0];
			sign_2 = dualgrid.vorticity_signs[3*vertical_index + 1];
			sign_3 = dualgrid.vorticity_signs[3*vertical_index + 2];
			out_field[i] = (1/dualgrid.area[i])*(grid.normal_distance[index_1]*sign_1*in_field[index_1] + grid.normal_distance[index_2]*sign_2*in_field[index_2] + grid.normal_distance[index_3]*sign_3*in_field[index_3]);
		}
		else
		{
			int index_1, index_2, index_3, index_4, sign_1, sign_2, sign_3, sign_4;
			double distance_1, distance_2, distance_3, distance_4;
			sign_2 = dualgrid.h_curl_signs[4*h_index + 1];
			sign_3 = dualgrid.h_curl_signs[4*h_index + 2];
			distance_2 = grid.normal_distance[index_2];
			distance_3 = grid.normal_distance[index_3];
			index_2 = dualgrid.h_curl_indices[4*h_index + 1];
			index_3 = dualgrid.h_curl_indices[4*h_index + 2];
			if(h_index >= NUMBER_OF_DUAL_VECTORS_H && h_index < NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_H)
			{
				index_1 = dualgrid.h_curl_indices[4*h_index + 0];
				index_4 = dualgrid.h_curl_indices[4*h_index + 3];
				sign_1= dualgrid.h_curl_signs[4*h_index + 0];
				sign_4 = dualgrid.h_curl_signs[4*h_index + 3];
				distance_1 = grid.normal_distance[index_1];
				distance_4 = grid.normal_distance[index_4];
			}
			else if (h_index < NUMBER_OF_DUAL_VECTORS_H)
			{
				index_1 = dualgrid.h_curl_indices[4*h_index + 3];
				index_4 = dualgrid.h_curl_indices[4*h_index + 3];
				sign_1 = -dualgrid.h_curl_signs[4*h_index + 3];
				sign_4 = dualgrid.h_curl_signs[4*h_index + 3];
				distance_1 = grid.normal_distance[index_1];
				distance_4 = grid.normal_distance[index_4];
			}
			else if (h_index >= NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_H)
			{
				index_1 = dualgrid.h_curl_indices[4*h_index + 0];
				index_4 = dualgrid.h_curl_indices[4*h_index + 0];
				sign_1 = dualgrid.h_curl_signs[4*h_index + 0];
				sign_4 = -dualgrid.h_curl_signs[4*h_index + 0];
				distance_1 = grid.normal_distance[index_1];
				distance_4 = grid.normal_distance[index_4];
				index_4 = dualgrid.h_curl_indices[4*h_index + 0];
			}
		out_field[i] = (1/dualgrid.area[i])*(distance_1*sign_1*in_field[index_1] + distance_2*sign_2*in_field[index_2] + distance_3*sign_3*in_field[index_3] + distance_4*sign_4*in_field[index_4]);
		}
	}
}
