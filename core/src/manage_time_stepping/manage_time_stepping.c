/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "../integrators/integrators.h"
#include "manage_time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int manage_time_stepping(State *state_old, State *state_new, Interpolate_info *interpolation, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, double delta_t)
{
	/*
	Here, the RK3 scheme is implemented.
	If radiation is updated, it is done at the first step.
	*/
	double delta_t_rk;
	for (int i = 0; i < 3; ++i)
	{
		// If tracers are on, phase transitions are only updated at the third step.
		config_info -> phase_transitions_on = 0;
		if (i == 2 && config_info -> tracers_on == 1)
		{
			config_info -> phase_transitions_on = 1;
		}
		delta_t_rk = delta_t/(3 - i);
		// Calculatung explicit horizontal momentum tendencies and vertical advective momentum tendencies.
		if (i == 0)
		{
			forward_tendencies(state_old, state_tendency, grid, dualgrid, diagnostics, forcings, interpolation, diffusion_info, config_info, i);
		}
		else
		{
			forward_tendencies(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, interpolation, diffusion_info, config_info, i);
		}
		// calculating the new values of the horizontal momentum, vertical horizontal advection handled implcitly
		three_band_solver_hor_vel_adv(state_old, state_new, state_tendency, delta_t_rk, grid);
		// Horizontal velocities can be considered as updated from now on.
		// The advective part of the vertical velocity equation is solved here, the vertical advection is calculated implicitly.
		three_band_solver_ver_vel_adv(state_old, state_new, state_tendency, delta_t_rk, grid);
		// here, the horizontal divergences are calculated with the new values of the horizontal velocity
		backward_tendencies(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t_rk, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, i);
		// determining the explicit component of the new temperature
		temperature_diagnostics_explicit(state_old, state_tendency, diagnostics, delta_t_rk);
		// here, the non-advective part of the vertical velocity equation is solved implicitly (sound wave solver)
		three_band_solver_ver_sound_waves(state_old, state_new, state_tendency, diagnostics, delta_t_rk, grid);
		// Vertical velocities can be seen as updated from now on.
		// now that the new vertical velocity is known, the new dry density can be calculated via implicit vertical advection
		three_band_solver_ver_den_dry(state_old, state_new, state_tendency, delta_t_rk, grid);
		// Now the entropy density is at the new step is calculated using the advection equation in flux form.
		three_band_solver_ver_entropy_density_dry(state_old, state_new, state_tendency, delta_t_rk, grid);
		// Vertical tracer advection with 3-band matrices.
		if (config_info -> tracers_on == 1)
		{
			three_band_solver_ver_tracers(state_old, state_new, state_tendency, delta_t_rk, grid);
			three_band_solver_ver_entropy_density_gaseous_tracers(state_old, state_new, state_tendency, delta_t_rk, grid);
		}
		// Diagnozing the temperature out of the linearized equation of state for energetic consistency.
		temperature_diagnostics(state_old, state_new);
    }
    return 0;
}








