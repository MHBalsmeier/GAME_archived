#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"

int is_vertical(int);

void rot(Vector_field in_field, Dual_vector_field out_field)
{
	extern Dualgrid dualgrid;
	extern Grid grid;
	for (int i = 0; i<=NUMBER_OF_DUAL_VECTORS-1; ++i)
	{
		if(dualgrid.is_vertical[i]>0)
		{
			int vertical_index;
			vertical_index = dualgrid.vertical_horizontal_index[i];
			int index_1, index_2, index_3, sign_1, sign_2, sign_3;
			index_1 = dualgrid.vorticity_indices[vertical_index][0];
			index_2 = dualgrid.vorticity_indices[vertical_index][1];
			index_3 = dualgrid.vorticity_indices[vertical_index][2];
			sign_1 = dualgrid.vorticity_signs[vertical_index][0];
			sign_2 = dualgrid.vorticity_signs[vertical_index][1];
			sign_3 = dualgrid.vorticity_signs[vertical_index][2];
			out_field[i] = (1/dualgrid.area[i])*(grid.normal_distance[index_1]*sign_1*in_field[index_1] +
			grid.normal_distance[index_2]*sign_2*in_field[index_2] + grid.normal_distance[index_3]*sign_3*in_field[index_3]);
		}
		else
		{
			int horizontal_index;
			horizontal_index = dualgrid.vertical_horizontal_index[i];
			int index_1, index_2, index_3, index_4, sign_1, sign_2, sign_3, sign_4;
			double distance_1, distance_2, distance_3, distance_4;
			sign_2 = dualgrid.h_curl_signs[horizontal_index][1];
			sign_3 = dualgrid.h_curl_signs[horizontal_index][2];
			distance_2 = grid.normal_distance[index_2];
			distance_3 = grid.normal_distance[index_3];
			index_2 = dualgrid.h_curl_indices[horizontal_index][1];
			index_3 = dualgrid.h_curl_indices[horizontal_index][2];
			if(horizontal_index >= NUMBER_OF_DUAL_VECTORS_H && horizontal_index < NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_H)
			{
				index_1 = dualgrid.h_curl_indices[horizontal_index][0];
				index_4 = dualgrid.h_curl_indices[horizontal_index][3];
				sign_1= dualgrid.h_curl_signs[horizontal_index][0];
				sign_4 = dualgrid.h_curl_signs[horizontal_index][3];
				distance_1 = grid.normal_distance[index_1];
				distance_4 = grid.normal_distance[index_4];
			}
			else if (horizontal_index < NUMBER_OF_DUAL_VECTORS_H)
			{
				index_1 = dualgrid.h_curl_indices[horizontal_index][3];
				index_4 = dualgrid.h_curl_indices[horizontal_index][3];
				sign_1 = -dualgrid.h_curl_signs[horizontal_index][3];
				sign_4 = dualgrid.h_curl_signs[horizontal_index][3];
				distance_1 = grid.normal_distance[index_1];
				distance_4 = grid.normal_distance[index_4];
			}
			else if (horizontal_index >= NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_H)
			{
				index_1 = dualgrid.h_curl_indices[horizontal_index][0];
				index_4 = dualgrid.h_curl_indices[horizontal_index][0];
				sign_1 = dualgrid.h_curl_signs[horizontal_index][0];
				sign_4 = -dualgrid.h_curl_signs[horizontal_index][0];
				distance_1 = grid.normal_distance[index_1];
				distance_4 = grid.normal_distance[index_4];
				index_4 = dualgrid.h_curl_indices[horizontal_index][0];
			}
		out_field[i] = (1/dualgrid.area[i])*(distance_1*sign_1*in_field[index_1] +
		distance_2*sign_2*in_field[index_2] + distance_3*sign_3*in_field[index_3] + 
		distance_4*sign_4*in_field[index_4]);

		}
	}
}
