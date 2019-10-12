#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"

void scalar_product(Vector_field in_field_1, Vector_field in_field_2, Scalar_field out_field)
{
	extern Grid grid;
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		int index_1, index_2, index_3, index_4, index_5,index_6, index_7, index_8;
		index_1 = grid.adjacent_vector_indices_h[i][0];
		index_2 = grid.adjacent_vector_indices_h[i][1];
		index_3 = grid.adjacent_vector_indices_h[i][2];
		index_4 = grid.adjacent_vector_indices_h[i][3];
		index_5 = grid.adjacent_vector_indices_h[i][4];
		int horizontal_index;
		horizontal_index = i-(i/NUMBER_OF_SCALARS_H)*NUMBER_OF_SCALARS_H;
		if (horizontal_index>NUMBER_OF_PENTAGONS-1)
			index_6 = grid.adjacent_vector_indices_h[i][5];
		index_7 = grid.adjacent_vector_index_lower[i];
		index_8 = grid.adjacent_vector_index_upper[i];
		out_field[i] = in_field_1[index_1]*in_field_2[index_1]
		+ in_field_1[index_2]*in_field_2[index_2]
		+ in_field_1[index_3]*in_field_2[index_3] + in_field_1[index_4]*in_field_2[index_4]
		+ in_field_1[index_5]*in_field_2[index_5] + in_field_1[index_6]*in_field_2[index_6]
		+ in_field_1[index_7]*in_field_2[index_7] + in_field_1[index_8]*in_field_2[index_8];
		if (horizontal_index>NUMBER_OF_PENTAGONS-1)
			out_field[i] = out_field[i] + in_field_1[index_6]*in_field_2[index_6];
	}
}
