/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int inner_product(Vector_field in_field_0, Vector_field in_field_1, Scalar_field out_field, Grid *grid, int vert_switch)
{
    // This function computes the inner product of the two vector fields in_field_0 and in_field_1. This is needed for computing the dissipation due to momentum diffusion (friction).
    int layer_index, h_index, j, no_of_edges;
    #pragma omp parallel for private (j, layer_index, h_index, no_of_edges)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
	    layer_index = i/NO_OF_SCALARS_H;
	    h_index = i - layer_index*NO_OF_SCALARS_H;
	    out_field[i] = 0;
        no_of_edges = 6;
        if (h_index < NO_OF_PENTAGONS)
        {
        	no_of_edges = 5;
        }
	    for (j = 0; j < no_of_edges; ++j)
	    {
	        out_field[i] += grid -> inner_product_weights[8*i + j]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]*
			in_field_1[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
	    }
	    out_field[i] += vert_switch*grid -> inner_product_weights[8*i + 6]*in_field_0[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]*in_field_1[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
	    out_field[i] += vert_switch*grid -> inner_product_weights[8*i + 7]*in_field_0[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]*in_field_1[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
	}
    return 0;
}

int kinetic_energy(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	// This function computes the kinetic energy, if in_field is the wind field. Only this part is neeed for the 3D Lamb transformation.
	inner_product(in_field, in_field, out_field, grid, 1);
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	out_field[i] = 0.5*out_field[i];
    }
    return 0;
}



