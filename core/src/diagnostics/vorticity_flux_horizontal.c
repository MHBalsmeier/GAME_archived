/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
This implements the modified TRSK scheme by Gassmann (2018).
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int vorticity_flux_horizontal(Vector_field in_field_0, Curl_field in_field_1, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    // From_index comes before to_index as usual.
	if (grid -> from_index[h_index] < NO_OF_PENTAGONS)
	{
		for (int i = 0; i < 4; ++i)
		{
			*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index + i]];
		}
	}
	else
	{
		for (int i = 0; i < 5; ++i)
		{
			if (i == 2)
			{
	    		*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*0.5*(in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index + i]] + in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + h_index]);
	    	}
	    	else
	    	{
	    		*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index + i]];
	    	}
		}
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
		for (int i = 5; i < 10; ++i)
		{
			if (i == 7)
			{
	    		*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*0.5*(in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index + i]] + in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + h_index]);
	    	}
	    	else
	    	{
	    		*component += grid -> trsk_modified_weights[10*h_index + i]*in_field_0[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]]*in_field_1[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index + i]];
	    	}
		}
	}
    return 0;
}










