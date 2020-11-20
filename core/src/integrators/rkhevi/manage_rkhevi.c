/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../integrators.h"
#include "../../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int manage_rkhevi(State *state_old, State *state_new, Interpolation_info *interpolation_info, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, double delta_t, double time_coordinate)
{
	/*
	Here, the RK3 scheme is implemented.
	If radiation is updated, it is done at the first step.
	*/
	double delta_t_rk;
	for (int i = 0; i < 3; ++i)
	{
		/*
		general remarks:
		-------------------------------------------------------------------------------
		At i == 0, it is state_new == state_old.
		state_old remains unchanged the whole time.
		*/
		
		// 1.) setting the time step of the RK substep
		delta_t_rk = delta_t/(3 - i);
		
		// 2.) Explicit component of the momentum equation.
		// ----------------------------------------------------------------------------
		forward_tendencies(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, interpolation_info, irreversible_quantities, config_info, i);
		
		// 3.) Vertical advection of horizontal momentum.
		// ---------------------------------------------------------------------------
		three_band_solver_ver_hor_vel_adv(state_old, state_tendency, state_new, delta_t, grid);
		
		// Horizontal velocity can be considered to be updated from now on.
		
		// 4.) Explicit component of the generalized density equations.
		// ----------------------------------------------------------------------------
		backward_tendencies(state_new, interpolation_info, state_tendency, grid, dualgrid, delta_t_rk, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, i, time_coordinate);
		// determining the explicit component of the new temperature
		
		// 5.) A pre-conditioned new temperature field, only containing explicit entropy and mass density tendencies (including diabatic forcings).
		// ----------------------------------------------------------------------------
		temperature_diagnostics_explicit(state_old, state_tendency, diagnostics, delta_t_rk);

		// 6.) Vertical sound wave solver.
		// ----------------------------------------------------------------------------
		three_band_solver_ver_sound_waves(state_old, state_tendency, state_new, diagnostics, delta_t_rk, grid);
		// Vertical velocity can be seen as updated from now on.
		
		// 7.) Solving the implicit component of the generalized density equaitons.
		// ----------------------------------------------------------------------------
		three_band_solver_gen_densitites(state_old, state_new, state_tendency, diagnostics, delta_t_rk, grid);
    }
    return 0;
}








