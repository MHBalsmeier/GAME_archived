/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include <stdio.h>

int recov_hor_ver_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    double vertical_factor, large_delta_z, small_delta_z;
    int aim_vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_VECTORS_V + h_index;
    large_delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 0]] - grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 2]];
    small_delta_z = grid -> z_vector[aim_vector_index] - grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 2]];
	vertical_factor = small_delta_z/large_delta_z;
	if (vertical_factor <= 0 || vertical_factor >= 1)
		printf("Error in recov_hor_ver_pri, position 0.\n");
    *component = 0.5*vertical_factor*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 0]];
	vertical_factor = 1 - vertical_factor;
    *component += 0.5*vertical_factor*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 2]];
    large_delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 1]] - grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 3]];
	small_delta_z = grid -> z_vector[aim_vector_index] - grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 3]];
	vertical_factor = small_delta_z/large_delta_z;
	if (vertical_factor <= 0 || vertical_factor >= 1)
		printf("Error in recov_hor_ver_pri, position 1.\n");
    *component += 0.5*vertical_factor*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 1]];
	vertical_factor = 1 - vertical_factor;
    *component += 0.5*vertical_factor*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> recov_hor_ver_pri_index[4*h_index + 3]];
    return 0;
}
