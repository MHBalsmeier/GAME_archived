/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, divergences get computed.
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include <stdio.h>

int divv_h(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This computes the divergence of a horizontal vector field.
	*/
	
    int layer_index, h_index, no_of_edges;
    double contra_upper, contra_lower, comp_h, comp_v;
	#pragma omp parallel for private(layer_index, h_index, no_of_edges, contra_upper, contra_lower, comp_h, comp_v)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        no_of_edges = 6;
        if (h_index < NO_OF_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        comp_h = 0;
        for (int j = 0; j < no_of_edges; ++j)
        {
			comp_h
			+= in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
			*grid -> adjacent_signs_h[6*h_index + j]
			*grid -> area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
        }
        comp_v = 0;
        if (layer_index == NO_OF_LAYERS - grid -> no_of_oro_layers - 1)
        {
            vertical_contravariant_corr(in_field, layer_index + 1, h_index, grid, &contra_lower);
            comp_v = -contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        else if (layer_index == NO_OF_LAYERS - 1)
        {
			vertical_contravariant_corr(in_field, layer_index, h_index, grid, &contra_upper);
			comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        }
        else if (layer_index > NO_OF_LAYERS - grid -> no_of_oro_layers - 1)
        {
            vertical_contravariant_corr(in_field, layer_index, h_index, grid, &contra_upper);
            vertical_contravariant_corr(in_field, layer_index + 1, h_index, grid, &contra_lower);
            comp_v
            = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]
            - contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        out_field[i] = 1/grid -> volume[i]*(comp_h + comp_v);
    }
    return 0;
}

int divv_h_limited(Vector_field in_field, Scalar_field out_field, Grid *grid, Scalar_field current_value, double delta_t)
{
	/*
	This is a limited version of the horizontal divergence operator.
	*/
	
	divv_h(in_field, out_field, grid);
	
    int layer_index, h_index, no_of_edges;
	double added_divergence, added_mass_rate;
	#pragma omp parallel for private(layer_index, h_index, no_of_edges, added_divergence, added_mass_rate)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		if (current_value[i] - delta_t*out_field[i] < 0)
		{
		    layer_index = i/NO_OF_SCALARS_H;
		    h_index = i - layer_index*NO_OF_SCALARS_H;
		    no_of_edges = 6;
		    if (h_index < NO_OF_PENTAGONS)
		    {
		    	no_of_edges = 5;
		    }
			added_divergence = current_value[i]/delta_t - out_field[i];
			out_field[i] += added_divergence;
			added_mass_rate = -added_divergence*grid -> volume[i];
			for (int j = 0; j < no_of_edges; ++j)
			{
				;
			}
		}
	}
	return 0;
}

int add_vertical_divv(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This adds the divergence of the vertical component of a vector field to the input scalar field.	
	*/
	
    int layer_index, h_index;
    double contra_upper, contra_lower, comp_v;
	#pragma omp parallel for private (layer_index, h_index, contra_upper, contra_lower, comp_v)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
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




