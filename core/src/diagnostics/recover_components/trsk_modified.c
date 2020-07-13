/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at horizontal vector points. If layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS, a vertical interpolation of the input field onto the desired z coordinate must be done.
*/

#include "../../enum_and_typedefs.h"
#include <stdio.h>

int trsk_modified(Vector_field in_field_0, Curl_field in_field_1, int layer_index, int h_index, double *component, Grid *grid)
{
	double unmodified_value, delta_z, dinputdz, delta_z_for_gradient, d_input;
    *component = 0;
    if (layer_index < NO_OF_LAYERS - NO_OF_ORO_LAYERS)
    {
		for (int i = 0; i < 10; ++i)
		    *component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]];
    }
    else
    {
		for (int i = 0; i < 10; ++i)
		{
			// vertical interpolation necessary
		    unmodified_value = grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]];
			if (layer_index == NO_OF_LAYERS - 1)
			{
				d_input = in_field_0[NO_OF_VECTORS_V + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]] - in_field_0[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]];
				delta_z_for_gradient = grid -> z_vector[NO_OF_VECTORS_V + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]] - grid -> z_vector[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			}
			else
			{
				d_input = in_field_0[NO_OF_VECTORS_V + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]] - in_field_0[NO_OF_VECTORS_V + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + grid -> trsk_modified_curl_indices[10*h_index + i]];
				delta_z_for_gradient = grid -> z_vector[NO_OF_VECTORS_V + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]] - grid -> z_vector[NO_OF_VECTORS_V + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			}
			dinputdz = d_input/delta_z_for_gradient;
			delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_VECTORS_V + h_index] - grid -> z_vector[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
			// result
		    *component += unmodified_value + delta_z*dinputdz;
	    }
    }
    return 0;
}










