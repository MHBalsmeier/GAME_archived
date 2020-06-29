/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at vertical vector points. The y component of a Vector_field is diagnozed here. If layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS, a vertical interpolation of the input field onto the desired z coordinate must be done.
*/

#include "../../enum_and_typedefs.h"

int recov_ver_1_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
    double z_aim, z_upper, z_lower, z_upper_weight, original_value, modified_value, delta_z, vertical_gradient, upper_value, delta_z_for_gradient;
    if (layer_index < NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
    {
		if (layer_index == 0)
		{
		    for (int i = 0; i < 6; ++i)
		        *component += grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field[NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		}
		else
		{
		    for (int i = 0; i < 6; ++i)
		    {
		        *component += 0.5*grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		        *component += 0.5*grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		    }
		}
    }
    else
    {
		if (layer_index == NUMBER_OF_LAYERS)
		{
		    for (int i = 0; i < 6; ++i)
		    {
		    	original_value = in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		    	upper_value = in_field[(layer_index - 2)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		    	delta_z_for_gradient = grid -> z_vector[(layer_index - 2)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]] - grid -> z_vector[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		    	vertical_gradient = (upper_value - original_value)/delta_z_for_gradient;
		    	delta_z = grid -> z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + h_index] - grid -> z_vector[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		    	modified_value = original_value + delta_z*vertical_gradient;
		        *component += grid -> recov_ver_1_pri_weight[6*h_index + i]*modified_value;
	        }
		}
		else
		{
		    for (int i = 0; i < 6; ++i)
		    {
				// vertical interpolation necessary
				z_aim = grid -> z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + h_index];
				z_upper = grid -> z_vector[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
				z_lower = grid -> z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
				z_upper_weight = (z_aim - z_lower)/(z_upper - z_lower);
		        *component += z_upper_weight*grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field[(layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
		        *component += (1 - z_upper_weight)*grid -> recov_ver_1_pri_weight[6*h_index + i]*in_field[layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + grid -> recov_ver_index[6*h_index + i]];
	        }
		}
    }
    return 0;
}















