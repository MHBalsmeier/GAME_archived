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
#include "atmostracers.h"
#include "../../diagnostics/diagnostics.h"
#include "../integrators.h"

int forward_tendencies(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Extrapolation_info *extrapolation_info, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_step_rk)
{
	// Update of the pressure gradient.
	manage_pressure_gradient(state, grid, dualgrid, diagnostics, forcings, extrapolation_info, irreversible_quantities, config_info, no_step_rk);
    if (no_step_rk == 0 && config_info -> momentum_diff == 1)
    {
		momentum_diff_diss(state, diagnostics, irreversible_quantities, config_info, grid, dualgrid);
		// Due to condensates, the friction acceleration needs to get a deceleration factor.
		if (config_info -> assume_lte == 0)
		{
			scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, irreversible_quantities -> friction_acc, irreversible_quantities -> friction_acc, grid);
		}
    }
	// Only the horizontal momentum is a forward tendency.
	integrate_momentum(state, state_tendency, grid, dualgrid, diagnostics, forcings, irreversible_quantities, config_info);
    return 0;
}
















