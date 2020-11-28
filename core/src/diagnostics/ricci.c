/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "diagnostics.h"

int vertical_contravariant(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	// Attention: adjacent_signs_h appears twice, thus does not need to be taken into account.
	if (h_index < 0 || h_index >= NO_OF_SCALARS_H)
	{
		return 1;
	}
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
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
				*result += -0.5*grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*in_field[vector_index];
			}
    	}
    	else if (layer_index == NO_OF_LAYERS)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result += -grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*in_field[vector_index];
			}
    	}
    	else
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result += -0.5*grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*in_field[vector_index];
			}
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result += -0.5*grid -> inner_product_weights[8*scalar_index + i]*grid -> slope[vector_index]*in_field[vector_index];
			}
    	}
    }
	return 0;
}

int horizontal_covariant(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	double vertical_component;
	remap_verpri2horpri_vector(in_field, layer_index, h_index, &vertical_component, grid);
	int vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_SCALARS_H + h_index;
	*result = in_field[vector_index] + grid -> slope[vector_index]*vertical_component;
	return 0;
}













