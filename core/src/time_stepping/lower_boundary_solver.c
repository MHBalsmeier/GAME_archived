/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdio.h>
#include <stdlib.h>

int solve_lower_boundary(State *state, Grid *grid)
{
    int layer_index = NO_OF_LAYERS;
	int layer_index_oro = layer_index - (NO_OF_LAYERS - NO_OF_ORO_LAYERS);
    double u_lowest, v_lowest, n_x, n_y, n_z, check_value;
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		n_x = grid -> vertical_contravar_unit[3*(layer_index_oro*NO_OF_SCALARS_H + h_index) + 0];
		n_y = grid -> vertical_contravar_unit[3*(layer_index_oro*NO_OF_SCALARS_H + h_index) + 1];
		n_z = grid -> vertical_contravar_unit[3*(layer_index_oro*NO_OF_SCALARS_H + h_index) + 2];			
		recov_ver_0_pri(state -> velocity_gas, layer_index, h_index, &u_lowest, grid);
		recov_ver_1_pri(state -> velocity_gas, layer_index, h_index, &v_lowest, grid);
		state -> velocity_gas[layer_index*NO_OF_VECTORS_PER_LAYER + h_index] = -1/n_z*(n_x*u_lowest + n_y*v_lowest);
		vertical_contravariant_normalized(state -> velocity_gas, layer_index, h_index, grid, &check_value);
		if (fabs(check_value) > 0.001)
		{
			printf("Error with lower boundary condition.\n");
			exit(1);	
		}
	}
	return 0;
}
