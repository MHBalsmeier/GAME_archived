/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../diagnostics.h"
#include <stdio.h>
#include <stdlib.h>

int solve_lower_boundary(State *state, Grid *grid)
{
	// int layer_index_oro = NO_OF_LAYERS - (NO_OF_LAYERS - NO_OF_ORO_LAYERS);
	double check_value, result_of_wind_h;
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		vertical_contravariant_normalized_h(state -> velocity_gas, NO_OF_LAYERS, h_index, grid, &result_of_wind_h);
		state -> velocity_gas[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + h_index] = -result_of_wind_h;
		check_value = 0;
		if (fabs(check_value) > 0.001)
		{
			printf("Error with lower boundary condition.\n");
			exit(1);	
		}
	}
	return 0;
}
