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

int forward_tendencies(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolation_info *interpolation_info, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_step_rk)
{
	// The calculation of the pressure gradient is only done at the first RK step.
	manage_pressure_gradient(current_state, grid, dualgrid, diagnostics, forcings, interpolation_info, irreversible_quantities, config_info, no_step_rk);
    if (no_step_rk == 2 && config_info -> momentum_diff == 1)
    {
		momentum_diff_diss(current_state, diagnostics, irreversible_quantities, config_info, grid, dualgrid);
		// Due to condensates, the friction acceleration needs to get a deceleration factor.
		scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, irreversible_quantities -> friction_acc, irreversible_quantities -> friction_acc, grid);
    }
	// Only the horizontal momentum is a forward tendency.
	integrate_momentum(current_state, state_tendency, grid, dualgrid, diagnostics, forcings, irreversible_quantities, config_info, no_step_rk);
    return 0;
}
















