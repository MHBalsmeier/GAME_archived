/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "diagnostics.h"

int vertical_contravariant_corr(Vector_field vector_field, int layer_index, int h_index, Grid *grid, double *result)
{
	/*
	Calculates (the vertical contravariant component - the vertical covariant component)
	of a vector field out of the horizontal contravariant components.
	*/
	// Attention: adjacent_signs_h appears twice, thus does not need to be taken into account.
	*result = 0;
	int scalar_index, vector_index;
	int no_of_edges = 6;
	double upper_volume, lower_volume, total_volume, vert_weight;
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
				lower_volume = grid -> volume_ratios[2*scalar_index + 0]*grid -> volume[scalar_index];
				upper_volume = grid -> volume_ratios[2*(scalar_index - NO_OF_SCALARS_H) + 1]*grid -> volume[scalar_index - NO_OF_SCALARS_H];
				total_volume = lower_volume + upper_volume;
				vert_weight = lower_volume/total_volume;
				*result += -vert_weight*grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*vector_field[vector_index];
			}
    	}
    	else if (layer_index == NO_OF_LAYERS)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result += -grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*vector_field[vector_index];
			}
    	}
    	else
    	{
			scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
			lower_volume = grid -> volume_ratios[2*scalar_index + 0]*grid -> volume[scalar_index];
			upper_volume = grid -> volume_ratios[2*(scalar_index - NO_OF_SCALARS_H) + 1]*grid -> volume[scalar_index - NO_OF_SCALARS_H];
			total_volume = lower_volume + upper_volume;
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				vert_weight = upper_volume/total_volume;
				*result += -vert_weight*grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*vector_field[vector_index];
			}
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				vert_weight = lower_volume/total_volume;
				*result += -vert_weight*grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*vector_field[vector_index];
			}
    	}
    }
	return 0;
}

int horizontal_covariant(Vector_field vector_field, int layer_index, int h_index, Grid *grid, double *result)
{
	// Calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components.
	double vertical_component;
	remap_verpri2horpri_vector(vector_field, layer_index, h_index, &vertical_component, grid);
	int vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
	*result = vector_field[vector_index] + grid -> slope[vector_index]*vertical_component;
	return 0;
}













