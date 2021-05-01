/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This file contains functions that perform averagings.
*/

#include "../enum_and_typedefs.h"
#include "geos95.h"
#include <stdio.h>

int remap_verpri2horpri_vector(Vector_field vector_field, int layer_index, int h_index, double *component, Grid *grid)
{
	/*
	reconstructs the vertical vector component *component at edge h_index in layer layer_index
	*/
    *component
    // layer above
    = grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + 6]
    *vector_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    *component
    += grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]) + 6]
    *vector_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
	// layer below
    if (layer_index < NO_OF_LAYERS - 1)
    {
		*component
		+= grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + 7]
		*vector_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
		*component
		+= grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]) + 7]
		*vector_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
    }
	// horizontal average
    *component = 0.5*(*component);
    return 0;
}

int vertical_contravariant_corr(Vector_field vector_field, int layer_index, int h_index, Grid *grid, double *result)
{
	/*
	calculates (the vertical contravariant component - the vertical covariant component)
	of a vector field out of the horizontal contravariant components
	*/
	// Attention: adjacent_signs_h appears twice, thus does not need to be taken into account.
	*result = 0;
	int scalar_index, vector_index;
	int no_of_edges = 6;
	if (h_index < NO_OF_PENTAGONS)
	{
		no_of_edges = 5;
	}
    if (layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
    {
    	if (layer_index == NO_OF_LAYERS - grid -> no_of_oro_layers)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -0.5
				*grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
    	}
    	else if (layer_index == NO_OF_LAYERS)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
    	}
    	else
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -0.5
				*grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -0.5
				*grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
    	}
    }
	return 0;
}

int horizontal_covariant(Vector_field vector_field, int layer_index, int h_index, Grid *grid, double *result)
{
	/*
	calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
	*/
	double vertical_component = 0;
	if (layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
	{
		remap_verpri2horpri_vector(vector_field, layer_index, h_index, &vertical_component, grid);
	}
	int vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
	*result = vector_field[vector_index] + grid -> slope[vector_index]*vertical_component;
	return 0;
}

int vector_field_hor_cov_to_con(Vector_field cov_to_con_field, Grid *grid)
{
    /*
    This function transforms the covariant horizontal measure numbers of a vector field to
    contravariant measure numbers.
    */
    int layer_index, h_index;
    double vertical_gradient;
	#pragma omp parallel for private(layer_index, h_index, vertical_gradient)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H && layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
        {
        	remap_verpri2horpri_vector(cov_to_con_field, layer_index, h_index - NO_OF_SCALARS_H, &vertical_gradient, grid);
            cov_to_con_field[i] += -grid -> slope[i]*vertical_gradient;
        }
    }
    return 0;
}

int tangential_wind(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
	/*
	This function computes the tangential component *component of the vector field in_field at edge h_index in layer layer_index
	using the TRSK weights.
	*/
	// initializing the result with zero
    *component = 0;
    // loop over the maximum of ten edges 
	for (int i = 0; i < 10; ++i)
	{
		*component += grid -> trsk_weights[10*h_index + i]
		*in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index + i]];
	}
    return 0;
}

int calc_uv_at_edge(Vector_field in_field, Vector_field out_field_u, Vector_field out_field_v, Grid *grid)
{
	/*
	This function diagnozes eastward and northward components of a vector field at edges.
	*/
	int layer_index, h_index;
	// orthogonal and tangential component at edge, respectively
	double wind_0, wind_1;
	#pragma omp parallel for private(layer_index, h_index, wind_0, wind_1)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		wind_0 = in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
		// finding the tangential component
		tangential_wind(in_field, layer_index, h_index, &wind_1, grid);
		// turning the Cartesian coordinate system to obtain u and v
		passive_turn(wind_0, wind_1, -grid -> direction[h_index],
		&out_field_u[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index],
		&out_field_v[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index]);
    }
	return 0;
}

int curl_field_to_cells(Curl_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This function averages a curl field from edges to cell centers.
	*/
	int layer_index, h_index, j, no_of_edges;
	#pragma omp parallel for private (j, layer_index, h_index, no_of_edges)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	// initializing the result with zero
        out_field[i] = 0;
        // determining the number of edges of the cell at hand
        no_of_edges = 6;
        if (h_index < NO_OF_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        // loop over all edges of the respective cell
        for (j = 0; j < no_of_edges; ++j)
        {
        	out_field[i] += 0.5
        	*grid -> inner_product_weights[8*i + j]
        	*in_field[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
    	}
    }
    return 0;
}

int edges_to_cells(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This function averages a vector field from edges to cell centers.
	*/
	int layer_index, h_index, j, no_of_edges;
	#pragma omp parallel for private (j, layer_index, h_index, no_of_edges)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
        // initializing the result with zero
        out_field[i] = 0;
        // determining the number of edges of the cell at hand
        no_of_edges = 6;
        if (h_index < NO_OF_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        // loop over all cell edges
        for (j = 0; j < no_of_edges; ++j)
        {
        	out_field[i] += 0.5
        	*grid -> inner_product_weights[8*i + j]
        	*in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
    	}
    }
    return 0;
}



