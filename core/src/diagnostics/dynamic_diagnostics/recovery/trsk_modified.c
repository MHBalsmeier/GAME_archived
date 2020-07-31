/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at horizontal vector points. If layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS, a vertical interpolation of the input field onto the desired z coordinate must be done (not yet implemented).
*/

#include "../../../enum_and_typedefs.h"
#include <stdio.h>

int trsk_modified(Vector_field in_field_0, Curl_field in_field_1, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    // From_index comes before to_index as usual.
	if (grid -> from_index[h_index] < NO_OF_PENTAGONS)
	{
		for (int i = 0; i < 4; ++i)
		{
			*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + i]];
		}
	}
	else
	{
	    *component += 
		grid -> trsk_modified_weights[10*h_index + 0]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 0]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 0]]
		+ grid -> trsk_modified_weights[10*h_index + 1]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 1]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 1]]
		+ grid -> trsk_modified_weights[10*h_index + 2]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 2]]*0.5*(in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 2]] + in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + h_index])
		+ grid -> trsk_modified_weights[10*h_index + 3]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 3]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 3]]
		+ grid -> trsk_modified_weights[10*h_index + 4]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 4]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 4]];
	}
	if (grid -> to_index[h_index] < NO_OF_PENTAGONS)	
	{
		for (int i = 5; i < 9; ++i)
		{
			*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index + i]];
		}
	}
	else
	{
	    *component += 
		grid -> trsk_modified_weights[10*h_index + 5]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 5]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 5]]
		+ grid -> trsk_modified_weights[10*h_index + 6]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 6]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 6]]
		+ grid -> trsk_modified_weights[10*h_index + 7]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 7]]*0.5*(in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 7]] + in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + h_index])
		+ grid -> trsk_modified_weights[10*h_index + 8]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 8]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 8]]
		+ grid -> trsk_modified_weights[10*h_index + 9]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + 9]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H  + grid -> trsk_modified_curl_indices[10*h_index + 9]];
	}
    return 0;
}










