/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/
/*
In this source file, the forward part of the integration is managed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "atmostracers.h"
#include "../diagnostics/diagnostics.h"
#include "time_stepping.h"

int forward_tendencies(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Extrapolation_info *extrapolation_info, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_step_rk, int slow_update_bool, double delta_t)
{
	// Update of the pressure gradient.
	if (no_step_rk == 0)
	{
		manage_pressure_gradient(state, grid, dualgrid, diagnostics, forcings, extrapolation_info, irreversible_quantities, config_info);
	}
	// Only the horizontal momentum is a forward tendency.
	vector_tendencies_expl(state, state_tendency, grid, dualgrid, diagnostics, forcings, irreversible_quantities, config_info, slow_update_bool, no_step_rk, delta_t);
    return 0;
}
















