/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdio.h>

int divv_h(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
    int layer_index, h_index, no_of_edges, i, j;
    double contra_upper, contra_lower, comp_h, comp_v;
	#pragma omp parallel for private(layer_index, h_index, no_of_edges, i, j, contra_upper, contra_lower, comp_h, comp_v)
    for (i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        no_of_edges = 6;
        if (h_index < NO_OF_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        comp_h = 0;
        for (j = 0; j < no_of_edges; ++j)
        {
			comp_h
			+= in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
			*grid -> adjacent_signs_h[6*h_index + j]
			*grid -> area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
        }
        if (layer_index == 0)
        {
			comp_v = 0;
        }
        else if (layer_index == NO_OF_LAYERS - 1)
        {
			vertical_contravariant_corr(in_field, layer_index, h_index, grid, &contra_upper);
			comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        }
        else
        {
            contra_upper = 0;
            if (layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
            {
                vertical_contravariant_corr(in_field, layer_index, h_index, grid, &contra_upper);
            }
            contra_lower = 0;
            if (layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers - 1)
            {
                vertical_contravariant_corr(in_field, layer_index + 1, h_index, grid, &contra_lower);
            }
            comp_v
            = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]
            - contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        out_field[i] = 1/grid -> volume[i]*(comp_h + comp_v);
    }
    return 0;
}

int add_vertical_divv(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This adds the divergence of the vertical component of a vector field to the input scalar field.	
	*/
    int layer_index, h_index, i;
    double contra_upper, contra_lower, comp_v;
	#pragma omp parallel for private (layer_index, h_index, i, contra_upper, contra_lower, comp_v)
    for (i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        if (layer_index == 0)
        {
        	contra_upper = 0;
        	contra_lower = in_field[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        else if (layer_index == NO_OF_LAYERS - 1)
        {
            contra_upper = in_field[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
            contra_lower = 0;
        }
        else
        {
            contra_upper = in_field[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
            contra_lower = in_field[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
    	comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]
    	- contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        out_field[i] += 1/grid -> volume[i]*comp_v;
    }
    return 0;
}




