/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this source file, the forward part of the integration is managed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../../manage_time_stepping/manage_time_stepping.h"
#include "atmostracers.h"
#include "../../diagnostics/diagnostics.h"
#include "../integrators.h"

int forward_tendencies(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	// The calculation of the pressure gradient is only done at the first RK step.
	if (no_step_rk == 0)
	{
		manage_pressure_gradient(current_state, grid, dualgrid, diagnostics, forcings, interpolation, diffusion_info, config_info, no_step_rk);
	}
    if (no_step_rk == 2)
    {
	if (config_info -> momentum_diff == 1)
	{
		momentum_diff_diss(current_state, diagnostics, diffusion_info, config_info, grid);
		// In the presence of condensates, the friction acceleration needs to get a deceleration factor.
		if (config_info -> tracers_on == 1)
		{
			scalar_times_vector(diffusion_info -> pressure_gradient_decel_factor, diffusion_info -> friction_acc, diffusion_info -> friction_acc, grid, 0);
		}
	}
    }
	// Only the horizontal momentum is a forward tendency.
	integrate_momentum(current_state, state_tendency, grid, dualgrid, diagnostics, forcings, diffusion_info, config_info, no_step_rk);
    return 0;
}
















