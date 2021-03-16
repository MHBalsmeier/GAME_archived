/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "diagnostics.h"
#include <stdio.h>

int vorticity_flux_vertical(Vector_field flux_density, Curl_field pot_vort, int layer_index, int h_index, double *component, Grid *grid, Dualgrid *dualgrid)
{
	/*
	Determining the vertical acceleration due to the vorticity flux term.
	*/
	// determining the number of edges
	int number_of_edges = 6;
	if (h_index < NO_OF_PENTAGONS)
	{
		number_of_edges = 5;
	}
	// initializing the result with zero
	*component = 0;
	// determining the vertical interpolation weight
	double vert_weight = 0.5;
	if (layer_index == 0 || layer_index == NO_OF_LAYERS)
	{
		vert_weight = 1;
	}
	if (layer_index >= 1)
	{
		for (int i = 0; i < number_of_edges; ++i)
		{
			*component +=
			vert_weight
			*grid -> inner_product_weights[8*((layer_index - 1)*NO_OF_SCALARS_H + h_index) + i]
			*flux_density[NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i]]
			*pot_vort[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + i]];
		}
	}
	if (layer_index <= NO_OF_LAYERS - 1)
	{
		for (int i = 0; i < number_of_edges; ++i)
		{
			*component +=
			vert_weight
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + h_index) + i]
			*flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i]]
			*pot_vort[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + i]];
		}
	}
    return 0;
}










