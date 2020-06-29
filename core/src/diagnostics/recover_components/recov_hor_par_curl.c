/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at horizontal vector points. The horizontal component of a Curl_field is diagnozed here. If layer_index != NO_OF_LAYERS - NO_OF_ORO_LAYERS && layer_index != NO_OF_LAYERS - NO_OF_ORO_LAYERS - 1, a vertical interpolation of the input field onto the desired z coordinate must be done.
*/

#include "../../enum_and_typedefs.h"

int recov_hor_par_curl(Curl_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
	double z_aim, z_upper, z_lower, z_weight_offset;
    *component = 0;
    if (layer_index != NO_OF_LAYERS - NO_OF_ORO_LAYERS && layer_index != NO_OF_LAYERS - NO_OF_ORO_LAYERS - 1)
    {
		for (int i = 0; i < 2; ++i)
		    *component += grid -> recov_hor_par_curl_weight[2*h_index + i]*in_field[layer_index*(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_hor_par_curl_index[2*h_index + i]];
    }
    else
    {
    	// vertical interpolation necessary
    	z_aim = grid -> z_vector[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
    	z_upper = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_par_curl_index[2*h_index + 0]];
    	z_lower = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_par_curl_index[2*h_index + 1]];
    	z_weight_offset = (z_aim - 0.5*(z_upper + z_lower))/(z_upper - z_lower);
    	// result (first index upper, second index lower)
    	*component += (grid -> recov_hor_par_curl_weight[2*h_index + 0] + z_weight_offset)*in_field[layer_index*(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_hor_par_curl_index[2*h_index + 0]] + (grid -> recov_hor_par_curl_weight[2*h_index + 1] - z_weight_offset)*in_field[layer_index*(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H) + grid -> recov_hor_par_curl_index[2*h_index + 1]];
    }
    return 0;
}
