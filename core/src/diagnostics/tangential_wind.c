/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"

int tangential_wind(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
	for (int i = 0; i < 10; ++i)
	{
		*component += grid -> trsk_weights[10*h_index + i]*in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index + i]];
	}
    return 0;
}
