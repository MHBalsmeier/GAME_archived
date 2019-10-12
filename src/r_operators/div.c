#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"


void divergence(Vector_field in_field, Scalar_field out_field)
{
	extern Grid grid;
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		int horizontal_index;
		horizontal_index = i-(i/NUMBER_OF_SCALARS_H)*NUMBER_OF_SCALARS_H;
		if (horizontal_index>NUMBER_OF_PENTAGONS-1)
		{
			out_field[i] =
			- in_field[grid.adjacent_vector_indices_h[i][0]]*grid.adjacent_signs_h[i][0]*grid.area[grid.adjacent_vector_indices_h[i][0]]
			- in_field[grid.adjacent_vector_indices_h[i][1]]*grid.adjacent_signs_h[i][1]*grid.area[grid.adjacent_vector_indices_h[i][1]]
			- in_field[grid.adjacent_vector_indices_h[i][2]]*grid.adjacent_signs_h[i][2]*grid.area[grid.adjacent_vector_indices_h[i][2]]
			- in_field[grid.adjacent_vector_indices_h[i][3]]*grid.adjacent_signs_h[i][3]*grid.area[grid.adjacent_vector_indices_h[i][3]]
			- in_field[grid.adjacent_vector_indices_h[i][4]]*grid.adjacent_signs_h[i][4]*grid.area[grid.adjacent_vector_indices_h[i][4]]
			- in_field[grid.adjacent_vector_indices_h[i][5]]*grid.adjacent_signs_h[i][5]*grid.area[grid.adjacent_vector_indices_h[i][5]]
			- in_field[grid.adjacent_vector_index_lower[i]]*grid.area[grid.adjacent_vector_index_lower[i]];
			+ in_field[grid.adjacent_vector_index_upper[i]]*grid.area[grid.adjacent_vector_index_upper[i]];
		}
		else
		{
			out_field[i] =
			- in_field[grid.adjacent_vector_indices_h[i][0]]*grid.adjacent_signs_h[i][0]*grid.area[grid.adjacent_vector_indices_h[i][0]]
			- in_field[grid.adjacent_vector_indices_h[i][1]]*grid.adjacent_signs_h[i][1]*grid.area[grid.adjacent_vector_indices_h[i][1]]
			- in_field[grid.adjacent_vector_indices_h[i][2]]*grid.adjacent_signs_h[i][2]*grid.area[grid.adjacent_vector_indices_h[i][2]]
			- in_field[grid.adjacent_vector_indices_h[i][3]]*grid.adjacent_signs_h[i][3]*grid.area[grid.adjacent_vector_indices_h[i][3]]
			- in_field[grid.adjacent_vector_indices_h[i][4]]*grid.adjacent_signs_h[i][4]*grid.area[grid.adjacent_vector_indices_h[i][4]]
			- in_field[grid.adjacent_vector_index_lower[i]]*grid.area[grid.adjacent_vector_index_lower[i]];
			+ in_field[grid.adjacent_vector_index_upper[i]]*grid.area[grid.adjacent_vector_index_upper[i]];
		}
		out_field[i] = (1/grid.volume[i])*out_field[i];
	}
}
