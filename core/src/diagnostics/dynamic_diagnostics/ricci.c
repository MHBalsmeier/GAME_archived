/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include <stdio.h>
#include "../diagnostics.h"

int vertical_contravariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	if (h_index < 0 || h_index >= NO_OF_SCALARS_H)
		return 1;
	int vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
	*result = in_field[vector_index];
	double horizontal_component;
	vertical_contravariant_normalized_h(in_field, layer_index, h_index, grid, &horizontal_component);
	*result += horizontal_component;
	return 0;
}

int vertical_contravariant_normalized_h(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	if (h_index < 0 || h_index >= NO_OF_SCALARS_H)
		return 1;
	*result = 0;
	int vector_index;
	int no_of_edges = 6;
	double vertical_coord_gradient;
	if (h_index < NO_OF_PENTAGONS)
	{
		no_of_edges = 5;
	}
    if (layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
    {
    	if (layer_index == NO_OF_LAYERS - NO_OF_ORO_LAYERS)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i];
				vertical_coord_gradient = (grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> to_index[grid -> recov_ver_index[6*h_index + i]]] - grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> from_index[grid -> recov_ver_index[6*h_index + i]]])/grid -> normal_distance[vector_index];
				*result += -0.5*fabs(grid -> recov_ver_weight[6*h_index + i])*vertical_coord_gradient*in_field[vector_index];
			}
    	}
    	else if (layer_index == NO_OF_LAYERS)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i];
				vertical_coord_gradient = (grid -> z_scalar[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[grid -> recov_ver_index[6*h_index + i]]] - grid -> z_scalar[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[grid -> recov_ver_index[6*h_index + i]]])/grid -> normal_distance[vector_index];
				*result += -fabs(grid -> recov_ver_weight[6*h_index + i])*vertical_coord_gradient*in_field[vector_index];
			}
    	}
    	else
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i];
				vertical_coord_gradient = (grid -> z_scalar[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[grid -> recov_ver_index[6*h_index + i]]] - grid -> z_scalar[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[grid -> recov_ver_index[6*h_index + i]]])/grid -> normal_distance[vector_index];
				*result += -0.5*fabs(grid -> recov_ver_weight[6*h_index + i])*vertical_coord_gradient*in_field[vector_index];
			}
			for (int i = 0; i < no_of_edges; ++i)
			{
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i];
				vertical_coord_gradient = (grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> to_index[grid -> recov_ver_index[6*h_index + i]]] - grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> from_index[grid -> recov_ver_index[6*h_index + i]]])/grid -> normal_distance[vector_index];
				*result += -0.5*fabs(grid -> recov_ver_weight[6*h_index + i])*vertical_coord_gradient*in_field[vector_index];
			}
    	}
    }
	return 0;
}

int horizontal_covariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	double vertical_component;
	recov_hor_ver_pri(in_field, layer_index, h_index, &vertical_component, grid);
	int vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_SCALARS_H + h_index;
	*result = in_field[vector_index] + grid -> tangential_coord_gradient[vector_index]*vertical_component;
	return 0;
}













