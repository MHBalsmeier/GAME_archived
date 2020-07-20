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

int manage_time_stepping(State *state_0, State *state_p1, Interpolate_info *interpolation, double delta_t, Grid *grid, Dualgrid *dualgrid, int momentum_diffusion_on, int rad_update, int tracers_on, int diffusion_on, Scalar_field radiation_tendency, State *state_tendency, Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divv, Scalar_field temperature, Vector_field t_tilde_flux_density, Scalar_field t_tilde_flux_density_divv, Vector_field temp_gradient, Scalar_field specific_entropy, Curl_field pot_vort, Vector_field pressure_gradient_acc, Vector_field pot_vort_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field e_kin_h, Vector_field pressure_gradient_acc_1, Vector_field temperature_flux_density, Vector_field temp_gradient_times_c_h_p, Vector_field pressure_gradient_acc_old, Vector_field e_kin_h_grad, Scalar_field wind_field_divv_h, int totally_first_step_bool, Scalar_field entropy_gas_flux_density_divv, Vector_field entropy_gas_flux_density, Diffusion_info *diffusion_info)
{
	/*
	Here, the RK3 scheme is implemented.
	If radiation is updated, it is done at the first step.
	*/
	int phase_transitions_on = 0;
	double delta_t_rk;
	for (int i = 0; i < 3; ++i)
	{
		if (totally_first_step_bool == 1)
			setup_interpolation(state_0, interpolation);
		if (i == 2 && tracers_on == 1)
			phase_transitions_on = 1;
		delta_t_rk = delta_t/(3 - i);
		// calculatung explicit horizontal momentum tendencies and vertical advective momentum tendencies
		if (i == 0)
		{
			explicit_momentum_tendencies(state_0, state_tendency, grid, dualgrid, momentum_diffusion_on, tracers_on, temperature, temp_gradient, specific_entropy, pot_vort, pressure_gradient_acc, pot_vort_tend, specific_entropy_gradient, c_h_p_field, e_kin_h, pressure_gradient_acc_1, momentum_diffusion_on, temp_gradient_times_c_h_p, pressure_gradient_acc_old, i, e_kin_h_grad, mass_dry_flux_density, totally_first_step_bool, diffusion_info);
		}
		else
		{
			explicit_momentum_tendencies(state_p1, state_tendency, grid, dualgrid, momentum_diffusion_on, tracers_on, temperature, temp_gradient, specific_entropy, pot_vort, pressure_gradient_acc, pot_vort_tend, specific_entropy_gradient, c_h_p_field, e_kin_h, pressure_gradient_acc_1, momentum_diffusion_on, temp_gradient_times_c_h_p, pressure_gradient_acc_old, i, e_kin_h_grad, mass_dry_flux_density, totally_first_step_bool, diffusion_info);
		}
		// calculating the new values of the horizontal momentum, vertical horizontal advection handled implcitly
		three_band_solver_hor(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// horizontal velocities can be considered updated from now on
		// now that the new values of the horizontal velocity values are knoown, the lower boundary condition can be solved
		solve_lower_boundary(state_p1, grid);
		// the advective part of the vertical velocity equation is solved here, the vertical advection is calculated implicitly
		three_band_solver_ver_vel_adv(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// vertical velocities can be updated from now on
		// here, the horizontal divergences are calculated with the new values of the horizontal velocity
		calc_partially_implicit_divvs(state_0, state_p1, interpolation, state_tendency, grid, dualgrid, momentum_diffusion_on, momentum_diffusion_on, tracers_on, delta_t, diffusion_on, radiation_tendency, phase_transitions_on, mass_dry_flux_density, mass_dry_flux_density_divv, temperature, t_tilde_flux_density, t_tilde_flux_density_divv, temp_gradient, temperature_flux_density, diffusion_on, wind_field_divv_h, i, entropy_gas_flux_density_divv, entropy_gas_flux_density, diffusion_info);
		// here, the non-advective part of the vertical velocity equation is solved implicitly (sound wave solver)
		three_band_solver_ver_sound_waves(state_0, state_p1, state_tendency, pressure_gradient_acc_1, t_tilde_flux_density_divv, wind_field_divv_h, delta_t_rk, grid);
		// now that the new vertical velocity is known, the new dry density can be calculated via implicit vertical advection
		three_band_solver_ver_den_dry(state_0, state_p1, state_tendency, delta_t_rk, grid);
		three_band_solver_ver_entropy_gas(state_0, state_p1, state_tendency, delta_t_rk, grid);
		// exit(1);
		// vertical tracer advection with 3-band matrices
		if (tracers_on == 1)
			three_band_solver_ver_tracers(state_0, state_p1, state_tendency, delta_t_rk, grid);
    }
	update_interpolation(state_0, state_p1, interpolation);
    return 0;
}







