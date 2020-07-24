/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "../diagnostics/diagnostics.h"
#include <omp.h>
#include "../spatial_operators/spatial_operators.h"

int set_interpolated_temperature(State *state_old, State *state_new, Interpolate_info *interpolate, int totally_first_bool)
{
	double old_weight = R_D/C_D_P - 0.5;
	double new_weight = 1 - old_weight;
	if (totally_first_bool == 1)
	{
		old_weight = 1;
		new_weight = 0;
	}
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		interpolate -> temp_interpolate[i] = old_weight*state_old -> temp_gas[i] + new_weight*state_new -> temp_gas[i];
	}
	return 0;
}
