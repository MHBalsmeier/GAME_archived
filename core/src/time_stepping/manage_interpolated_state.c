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

int setup_interpolated_state(State *state_old, State_interpolate *state_interpolate)
{
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		state_interpolate -> t_tilde[i] = state_old -> t_tilde[i];
	}
	return 0;
}

int update_interpolated_state(State *state_old, State *state_new, State_interpolate *state_interpolate)
{
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		state_interpolate -> t_tilde[i] = state_new -> t_tilde[i] + 0.5*(state_new -> t_tilde[i] - state_old -> t_tilde[i]);
	}
	return 0;
}
