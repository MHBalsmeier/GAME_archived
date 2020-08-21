/*
This source file is part of the General Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
This function revocers vertical components of primal vectors at edges. Terrain following ccordinates are taken into account through vertical interpolation.
*/

#include "../../../enum_and_typedefs.h"
#include <stdio.h>

int recov_hor_ver_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0.5*grid -> volume_ratios[2*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + 0]*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    *component += 0.5*grid -> volume_ratios[2*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + 1]*in_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    *component += 0.5*grid -> volume_ratios[2*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]) + 0]*in_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
    *component += 0.5*grid -> volume_ratios[2*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]) + 1]*in_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
    return 0;
}
