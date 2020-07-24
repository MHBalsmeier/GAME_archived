/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at vertical vector points. The x component of a Curl_field is diagnozed here. If layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS, a vertical interpolation of the input field onto the desired z coordinate must be done.
*/

#include "../../../enum_and_typedefs.h"

int recov_ver_0_curl(Curl_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
	double unmodified_value, d_input, delta_z_for_gradient, delta_z, dinputdz;
    *component = 0;
    if (layer_index < NO_OF_LAYERS - NO_OF_ORO_LAYERS)
    {
		for (int i = 0; i < 6; ++i)
		    *component += grid -> recov_ver_0_curl_weight[6*h_index + i]*in_field[layer_index*(NO_OF_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]];
	}
	else
	{
		for (int i = 0; i < 6; ++i)
		{
			// vertical interpolation necessary
			unmodified_value = grid -> recov_ver_0_curl_weight[6*h_index + i]*in_field[layer_index*(NO_OF_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]];
			if (layer_index == NO_OF_LAYERS)
			{
				d_input = in_field[(layer_index - 1)*(NO_OF_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]] - in_field[layer_index*(NO_OF_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]];
				delta_z_for_gradient = grid -> z_vector[(layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i]] - grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i]];
			}
			else
			{
				d_input = in_field[(layer_index - 1)*(NO_OF_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]] - in_field[(layer_index + 1)*(NO_OF_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_ver_index[6*h_index + i]];
				delta_z_for_gradient = grid -> z_vector[(layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i]] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i]];
			}
			dinputdz = d_input/delta_z_for_gradient;
			delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + h_index] - grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_ver_index[6*h_index + i]];
			// result
		    *component += unmodified_value + delta_z*dinputdz;
	    }
	}
    return 0;
}
