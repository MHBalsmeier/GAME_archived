/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int setup_interpolation(State *state_old, Interpolate_info *interpolation)
{
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		interpolation -> t_tilde[i] = state_old -> t_tilde[i];
	}
	return 0;
}

int update_interpolation(State *state_old, State *state_new, Interpolate_info *interpolation)
{
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		interpolation -> t_tilde[i] = state_new -> t_tilde[i] + 0.5*(state_new -> t_tilde[i] - state_old -> t_tilde[i]);
	}
	return 0;
}
