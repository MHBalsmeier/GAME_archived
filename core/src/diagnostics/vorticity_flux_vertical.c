/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "diagnostics.h"
#include <stdio.h>

int vorticity_flux_vertical(Vector_field in_field_0, Curl_field in_field_1, int layer_index, int h_index, double *component, Grid *grid, Dualgrid *dualgrid)
{
	int number_of_edges = 6;
	if (h_index < NO_OF_PENTAGONS)
	{
		number_of_edges = 5;
	}
	double vector_field_value;
	*component = 0;
	int layer_index_mod = layer_index;
	if (layer_index == NO_OF_LAYERS)
	{
		layer_index_mod = layer_index - 1;
	}
	for (int i = 0; i < number_of_edges; ++i)
	{
		remap_horpri2hordual_vector(in_field_0, layer_index, grid -> adjacent_vector_indices_h[6*h_index + i], &vector_field_value, grid);
		*component +=
		// the inner product weights are used here to remap the quantities to the scalar data points
		grid -> inner_product_weights[8*(layer_index_mod*NO_OF_SCALARS_H + h_index) + i]
		*vector_field_value
		*in_field_1[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + i]];
	}
    return 0;
}










