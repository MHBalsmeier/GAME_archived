/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
This function revocers vertical components of primal vectors at edges. Terrain following ccordinates are taken into account through vertical interpolation.
*/

#include "../../../enum_and_typedefs.h"
#include <stdio.h>

int recov_primal2dual(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
	if (layer_index == 0)
	{
		*component = in_field[NO_OF_SCALARS_H + h_index];
	}
	else if (layer_index == NO_OF_LAYERS)
	{
		*component = in_field[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + h_index];
	}
	else
	{
		*component = grid -> recov_primal2dual_weights[2*(layer_index*NO_OF_VECTORS_H + h_index) + 0]*in_field[NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + h_index];
		*component += grid -> recov_primal2dual_weights[2*(layer_index*NO_OF_VECTORS_H + h_index) + 1]*in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
	}
    return 0;
}
