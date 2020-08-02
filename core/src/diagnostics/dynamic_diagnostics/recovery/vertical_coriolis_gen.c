/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../../enum_and_typedefs.h"
#include "../../diagnostics.h"
#include <stdio.h>

int vertical_coriolis_gen(Vector_field in_field_0, Curl_field in_field_1, int layer_index, int h_index, double *component, Grid *grid)
{
	int number_of_edges = 6;
	if (h_index < NO_OF_PENTAGONS)
		number_of_edges = 5;
	double vector_field_value;
	*component = 0;
	for (int i = 0; i < number_of_edges; ++i)
	{
		recov_primal2dual(in_field_0, layer_index, grid -> adjacent_vector_indices_h[6*h_index + i], &vector_field_value, grid);
		*component += grid -> recov_ver_weight[6*(layer_index*NO_OF_SCALARS_H + h_index) + i]*vector_field_value*in_field_1[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + i]];
	}
    return 0;
}










