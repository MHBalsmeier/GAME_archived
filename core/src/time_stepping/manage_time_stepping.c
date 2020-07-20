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

int manage_time_stepping(State *state_0, State *state_p1, Interpolate_info *interpolation, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, double delta_t)
{
	/*
	Here, the RK3 scheme is implemented.
	If radiation is updated, it is done at the first step.
	*/
	double delta_t_rk;
	for (int i = 0; i < 3; ++i)
	{
		if (config_info -> totally_first_step_bool == 1)
			setup_interpolation(state_0, interpolation);
		if (i == 2 && config_info -> tracers_on == 1)
			config_info -> phase_transitions_on = 1;
		delta_t_rk = delta_t/(3 - i);
		// Calculatung explicit horizontal momentum tendencies and vertical advective momentum tendencies.
		if (i == 0)
			explicit_momentum_tendencies(state_0, state_tendency, grid, dualgrid, diagnostics, forcings, interpolation, diffusion_info, config_info, i);
		else
			explicit_momentum_tendencies(state_p1, state_tendency, grid, dualgrid, diagnostics, forcings, interpolation, diffusion_info, config_info, i);
		// calculating the new values of the horizontal momentum, vertical horizontal advection handled implcitly
		three_band_solver_hor(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// Horizontal velocities can be considered as updated from now on.
		// Now that the new values of the horizontal velocity values are known, the lower boundary condition can be solved.
		solve_lower_boundary(state_p1, grid);
		// The advective part of the vertical velocity equation is solved here, the vertical advection is calculated implicitly.
		three_band_solver_ver_vel_adv(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// Vertical velocities can be seen as updated from now on.
		// here, the horizontal divergences are calculated with the new values of the horizontal velocity
		calc_partially_implicit_divvs(state_0, state_p1, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, i);
		// here, the non-advective part of the vertical velocity equation is solved implicitly (sound wave solver)
		three_band_solver_ver_sound_waves(state_0, state_p1, state_tendency, diagnostics, forcings, delta_t_rk, grid);
		// now that the new vertical velocity is known, the new dry density can be calculated via implicit vertical advection
		three_band_solver_ver_den_dry(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// Now the entropy density is at the new step is calculated using the advection equation in flux form.
		three_band_solver_ver_entropy_gas(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// exit(1);
		// Vertical tracer advection with 3-band matrices.
		if (config_info -> tracers_on == 1)
			three_band_solver_ver_tracers(state_0, state_p1, state_tendency, delta_t_rk, grid);
    }
	update_interpolation(state_0, state_p1, interpolation);
    return 0;
}







