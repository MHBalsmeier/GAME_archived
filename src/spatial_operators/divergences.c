/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, divergences get computed.
*/

#include <stdio.h>
#include "../game_types.h"
#include "spatial_operators.h"

int divv_h(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This computes the divergence of a horizontal vector field.
	*/
	
    int i, no_of_edges;
    double contra_upper, contra_lower, comp_h, comp_v;
	#pragma omp parallel for private(i, no_of_edges, contra_upper, contra_lower, comp_h, comp_v)
    for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
    {
	    no_of_edges = 6;
	    if (h_index < NO_OF_PENTAGONS)
	    {
	    	no_of_edges = 5;
	    }
    	for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
    	{
		    i = layer_index*NO_OF_SCALARS_H + h_index;
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
    }
    return 0;
}

int divv_h_limited(Vector_field in_field, Scalar_field out_field, Grid *grid, Scalar_field current_value, double delta_t)
{
	/*
	This is a limited version of the horizontal divergence operator.
	*/
	
	// calling the normal horizontal divergence operator
	divv_h(in_field, out_field, grid);
	
    int i, no_of_edges, adjacent_scalar_index_h;
	double added_divergence, added_mass_rate, out_flow_rate_sum, outflow_rate, outflow_rate_factor;
	#pragma omp parallel for private(i, no_of_edges, added_divergence, added_mass_rate, out_flow_rate_sum, outflow_rate, outflow_rate_factor, adjacent_scalar_index_h)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		no_of_edges = 6;
		if (h_index < NO_OF_PENTAGONS)
		{
			no_of_edges = 5;
		}
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			// Negative mass densities are not possible. This is the case we want to look at.
			if (current_value[i] - delta_t*out_field[i] < 0.0)
			{
				// this is the excess divergence (negative)
				added_divergence = current_value[i]/delta_t - out_field[i];
				// we add the excess divergence so that the operator produces no negative mass densities
				out_field[i] += added_divergence;
				// this is the additional mass source rate that is being produced by the limiter and that violates the mass conservation
				added_mass_rate = -added_divergence*grid -> volume[i];
				// determining how much mass flows out of the grid box per time interval
				out_flow_rate_sum = 0.0;
				for (int j = 0; j < no_of_edges; ++j)
				{
					outflow_rate = in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
					*grid -> adjacent_signs_h[6*h_index + j]
					*grid -> area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
					// only positive values are counted here because we are only interestedin what flows out of the grid box
					if (outflow_rate > 0.0)
					{
						out_flow_rate_sum += outflow_rate;
					}
				}
				// now we want to reduce what flows out of the grid box
				// outflow_rate_factor*out_flow_rate_sum = added_mass_rate
				outflow_rate_factor = added_mass_rate/out_flow_rate_sum;
				for (int j = 0; j < no_of_edges; ++j)
				{
					// rescaling everything that flows out
					if (in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]*grid -> adjacent_signs_h[6*h_index + j] > 0)
					{
						adjacent_scalar_index_h = grid -> to_index[grid -> adjacent_vector_indices_h[6*h_index  + j]];
						if (adjacent_scalar_index_h == h_index)
						{
							adjacent_scalar_index_h = grid -> from_index[grid -> adjacent_vector_indices_h[6*h_index  + j]];
						}
						out_field[layer_index*NO_OF_SCALARS_H + adjacent_scalar_index_h] += outflow_rate_factor
						*in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
						*grid -> adjacent_signs_h[6*h_index + j]
						*grid -> area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
						/grid -> volume[layer_index*NO_OF_SCALARS_H + adjacent_scalar_index_h];
					}
				}
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
	
    int i;
    double contra_upper, contra_lower, comp_v;
	#pragma omp parallel for private (i, contra_upper, contra_lower, comp_v)
    for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
    {
    	for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
    	{
    		i = layer_index*NO_OF_SCALARS_H + h_index;
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
    }
    return 0;
}




