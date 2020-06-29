/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at horizontal vector points. If layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS, a vertical interpolation of the input field onto the desired z coordinate must be done.
*/

#include "../../enum_and_typedefs.h"

int recov_hor_par_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
	double unmodified_value, delta_z, dinputdz, delta_z_for_gradient, d_input;
    *component = 0;
    if (layer_index < NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
    {
		for (int i = 0; i < 10; ++i)
		    *component += grid -> trsk_modified_weights[10*h_index + i]*in_field[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
    }
    else
    {
		for (int i = 0; i < 10; ++i)
		{
			// vertical interpolation necessary
		    unmodified_value = grid -> trsk_modified_weights[10*h_index + i]*in_field[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			if (layer_index == NUMBER_OF_LAYERS - 1)
			{
				d_input = in_field[NUMBER_OF_VECTORS_V + (layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]] - in_field[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
				delta_z_for_gradient = grid -> z_vector[NUMBER_OF_VECTORS_V + (layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]] - grid -> z_vector[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			}
			else
			{
				d_input = in_field[NUMBER_OF_VECTORS_V + (layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]] - in_field[NUMBER_OF_VECTORS_V + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
				delta_z_for_gradient = grid -> z_vector[NUMBER_OF_VECTORS_V + (layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]] - grid -> z_vector[NUMBER_OF_VECTORS_V + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			}
			dinputdz = d_input/delta_z_for_gradient;
			delta_z = grid -> z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + h_index] - grid -> z_vector[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			// result
		    *component += unmodified_value + delta_z*dinputdz;
	    }
    }
    return 0;
}
