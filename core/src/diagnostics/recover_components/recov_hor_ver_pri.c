/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
This function revocers vertical components of primal vectors at edges. Terrain following ccordinates are taken into account through vertical interpolation.
*/

#include "../../enum_and_typedefs.h"
#include <stdio.h>

int recov_hor_ver_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    double vertical_factor, large_delta_z, small_delta_z;
    // This is the index of the edge we aim at.
    int aim_vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_VECTORS_V + h_index;
    // Now comes the first neighbouring cell.
    large_delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    small_delta_z = grid -> z_vector[aim_vector_index] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
	vertical_factor = small_delta_z/large_delta_z;
	// The resulting component is a sum of four constituents.
    *component = 0.5*vertical_factor*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
	vertical_factor = 1 - vertical_factor;
    *component += 0.5*vertical_factor*in_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    // Now comes the other neighbouring cell.
    large_delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
	small_delta_z = grid -> z_vector[aim_vector_index] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
	vertical_factor = small_delta_z/large_delta_z;
    *component += 0.5*vertical_factor*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
	vertical_factor = 1 - vertical_factor;
    *component += 0.5*vertical_factor*in_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
    return 0;
}
