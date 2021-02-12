/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "diagnostics.h"
#include <stdio.h>

int vorticity_flux_vertical(Vector_field flux_density, Curl_field pot_vort, int layer_index, int h_index, double *component, Grid *grid, Dualgrid *dualgrid)
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
		remap_horpri2hordual_vector(flux_density, layer_index, grid -> adjacent_vector_indices_h[6*h_index + i], &vector_field_value, grid);
		*component +=
		// the inner product weights are used here to remap the quantities to the scalar data points
		grid -> inner_product_weights[8*(layer_index_mod*NO_OF_SCALARS_H + h_index) + i]
		*vector_field_value
		*pot_vort[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + i]];
	}
    return 0;
}










