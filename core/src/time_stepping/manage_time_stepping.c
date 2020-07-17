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

int manage_time_stepping(State *state_0, State *state_p1, State *state_interpolate, double delta_t, Grid *grid, Dualgrid *dualgrid, int momentum_diffusion_on, int rad_update, int tracers_on, int diffusion_on, Scalar_field radiation_tendency, double tracer_mass_source_rates[], double tracer_heat_source_rates[], State *state_tendency, Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divv, Scalar_field temperature, Vector_field entropy_gas_flux_density, Scalar_field entropy_gas_flux_density_divv, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Vector_field friction_acc, Scalar_field heating_diss, Scalar_field specific_entropy, Curl_field pot_vort, Vector_field gradient_geopotential_energy, Vector_field pressure_gradient_acc, Vector_field pot_vort_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field macroscopic_energy, Scalar_field pressure_gradient_decel_factor, Vector_field pressure_gradient_acc_1, Scalar_field diffusion_coeff_numerical_h, Scalar_field diffusion_coeff_numerical_v, Vector_field mass_dry_diffusion_flux_density, Scalar_field mass_dry_diffusion_source_rate, Vector_field temperature_flux_density, Scalar_field tracer_density, Vector_field tracer_velocity, Vector_field tracer_flux_density, Scalar_field tracer_flux_density_divv, Scalar_field tracer_density_temperature, Vector_field tracer_temperature_flux_density, Scalar_field tracer_temperature_flux_density_divv, Vector_field temp_gradient_times_c_h_p, Vector_field pressure_gradient_acc_old, int first_step_bool, Vector_field e_kin_h_grad, Scalar_field temperature_density, Scalar_field temperature_flux_density_divv, Scalar_field wind_field_divv, int totally_first_step_bool)
{
	/*
	Here, the RK3 scheme is implemented.
	If radiation is updated, it is done at the first step.
	*/
    linear_combine_two_states(state_0, state_0, state_interpolate, 1, 0);
	int phase_transitions_on = 0;
	double delta_t_rk;
	for (int i = 0; i < 3; ++i)
	{
		if (i == 2 && tracers_on == 1)
			phase_transitions_on = 1;
		delta_t_rk = delta_t/(3 - i);
		// calculatung explicit horizontal momentum tendencies and vertical advective momentum tendencies
		explicit_momentum_tendencies(state_interpolate, state_tendency, grid, dualgrid, momentum_diffusion_on, tracers_on, temperature, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, pot_vort, gradient_geopotential_energy, pressure_gradient_acc, pot_vort_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, momentum_diffusion_on, temp_gradient_times_c_h_p, pressure_gradient_acc_old, i, e_kin_h_grad, mass_dry_flux_density, totally_first_step_bool);
		// calculating the new values of the horizontal momentum, vertical horizontal advection handled implcitly
		three_band_solver_hor(state_interpolate, state_p1, state_tendency, delta_t, grid);
		// horizontal velocities can be considered updated from now on
		// now that the new values of the horizontal velocity values are knoown, the lower boundary condition can be solved
		solve_lower_boundary(state_p1, grid);
		// the advective part of the vertical velocity equation is solved here, the vertical advection is calculated implicitly
		three_band_solver_ver_vel_adv(state_interpolate, state_p1, state_tendency, delta_t, grid);
		// vertical velocities can be updated from now on
		// here, the horizontal divergences are calculated with the new values of the horizontal velocity
		calc_partially_implicit_divvs(state_interpolate, state_p1, state_tendency, grid, dualgrid, momentum_diffusion_on, momentum_diffusion_on, tracers_on, delta_t, diffusion_on, radiation_tendency, phase_transitions_on, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divv, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divv, temp_diffusion_heating, temp_gradient, heating_diss, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divv, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divv, diffusion_on);
		// here, the non-advective part of the vertical velocity equation is solved implicitly (sound wave solver)
		three_band_solver_ver_sound_waves(state_interpolate, state_p1, state_tendency, pressure_gradient_acc_1, gradient_geopotential_energy, temperature, temperature_density, temperature_flux_density, temperature_flux_density_divv, wind_field_divv, delta_t, grid);
		// now that the new vertical velocity is known, the new dry density can be calculated via implicit vertical advection
		three_band_solver_ver_den_dry(state_interpolate, state_p1, state_tendency, delta_t, grid);
		// Here, the global entropy density is calculated via its generalized equation of state.
		three_band_solver_ver_entropy_gas(state_interpolate, state_p1, state_tendency, delta_t, grid);
		// vertical tracer advection with 3-band matrices
		if (tracers_on == 1)
			three_band_solver_ver_tracers(state_interpolate, state_p1, state_tendency, delta_t, grid);
		// diagnozing the tendency
		linear_combine_two_states(state_interpolate, state_p1, state_tendency, -1/delta_t, 1/delta_t);
		// determining the new interpolated state
		if (i < 2)
			linear_combine_two_states(state_0, state_tendency, state_interpolate, 1, delta_t_rk);
    }
	// the final step
    linear_combine_two_states(state_0, state_tendency, state_p1, 1, delta_t);
    return 0;
}







