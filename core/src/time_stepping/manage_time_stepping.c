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

int manage_time_stepping(State *state_0, State *state_p1, State *state_p2, State *state_tendency_precond, double delta_t, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_update, int tracers_on, int diffusion_on, Scalar_field radiation_tendency, double tracer_mass_source_rates[], double tracer_heat_source_rates[], State *state_tendency, Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divv, Scalar_field temperature, Vector_field entropy_gas_flux_density, Scalar_field entropy_gas_flux_density_divv, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Vector_field friction_acc, Scalar_field heating_diss, Scalar_field specific_entropy, Curl_field rel_curl, Curl_field abs_curl, Vector_field gradient_geopotential_energy, Vector_field pressure_gradient_acc, Vector_field abs_curl_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field macroscopic_energy, Scalar_field pressure_gradient_decel_factor, Vector_field pressure_gradient_acc_1, Scalar_field diffusion_coeff_numerical_h, Scalar_field diffusion_coeff_numerical_v, Vector_field mass_dry_diffusion_flux_density, Scalar_field mass_dry_diffusion_source_rate, Vector_field temperature_flux_density, Scalar_field tracer_density, Vector_field tracer_velocity, Vector_field tracer_flux_density, Scalar_field tracer_flux_density_divv, Scalar_field tracer_density_temperature, Vector_field tracer_temperature_flux_density, Scalar_field tracer_temperature_flux_density_divv, Vector_field temp_gradient_times_c_h_p, Vector_field pressure_gradient_acc_old, int first_step_bool, Vector_field e_kin_h_grad, Scalar_field temperature_density, Scalar_field temperature_flux_density_divv, Scalar_field wind_field_divv)
{
	/*
	Here, the predictor-corrector is implemented (it is currently rather an RK2 scheme).
	Diffusion and dissipation are only updated at the predictor step.
	If radiation is updated, it is done at the predictor step.
	*/
	// predictor step
	// calculatung explicit horizontal momentum tendencies and vertical horizontal advective momentum tendencies
	explicit_momentum_tendencies(state_0, state_tendency_precond, grid, dualgrid, dissipation_on, tracers_on, temperature, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, rel_curl, abs_curl, gradient_geopotential_energy, pressure_gradient_acc, abs_curl_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, dissipation_on, temp_gradient_times_c_h_p, pressure_gradient_acc_old, 0, first_step_bool, e_kin_h_grad);
	// calculating the new values of the horizontal momentum, vertical advection handled implcitly
	three_band_solver_hor(state_0, state_p1, state_tendency_precond, 0.5*delta_t, grid);
	// now that the new values of the horizontal velocity values are knoown, the lower boundary condition can be solved
	solve_lower_boundary(state_p1, grid);
	// the advective part of the vertical velocity equation is solved here, the vertical advection is calculated implicitly
	three_band_solver_ver_vel_adv(state_0, state_p1, state_tendency_precond, 0.5*delta_t, grid);
	// here, the horizontal divergences are calculated with the new values of the horizontal velocity
	calc_partially_implicit_divvs(state_p1, state_tendency_precond, grid, dualgrid, dissipation_on, rad_update, tracers_on, delta_t, diffusion_on, radiation_tendency, tracers_on, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divv, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divv, temp_diffusion_heating, temp_gradient, heating_diss, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divv, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divv, diffusion_on);
	// here, the non-advective part of the vertical velocity equation is solved implicitly (sound wave solver)
	three_band_solver_ver_sound_waves(state_0, state_p1, state_tendency_precond, pressure_gradient_acc_1, gradient_geopotential_energy, temperature, temperature_density, temperature_flux_density, temperature_flux_density_divv, wind_field_divv, 0.5*delta_t, grid);
	// now that the new vertical velocity is known, the new dry density can be calculated via implicit vertical advection
	three_band_solver_ver_den_dry(state_0, state_p1, state_tendency_precond, 0.5*delta_t, grid);
	// Here, the global entropy density is calculated via its generalized equation of state.
	three_band_solver_ver_entropy_gas(state_0, state_p1, state_tendency_precond, 0.5*delta_t, grid);
	// diagnozing the tendency
	linear_combine_two_states(state_0, state_p1, state_tendency_precond, -1/(0.5*delta_t), 1/(0.5*delta_t));
	// vertical tracer advection with 3-band matrices
	if (tracers_on == 1)
		three_band_solver_ver_tracers(state_0, state_p1, state_tendency_precond, 0.5*delta_t, grid);
	// corrector step
	// calculatung explicit horizontal momentum tendencies and vertical advective momentum tendencies
	explicit_momentum_tendencies(state_p1, state_tendency, grid, dualgrid, dissipation_on, tracers_on, temperature, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, rel_curl, abs_curl, gradient_geopotential_energy, pressure_gradient_acc, abs_curl_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, 0, temp_gradient_times_c_h_p, pressure_gradient_acc_old, 1, first_step_bool, e_kin_h_grad);
	// calculating the new values of the horizontal momentum, vertical horizontal advection handled implcitly
	three_band_solver_hor(state_p1, state_p2, state_tendency, delta_t, grid);
	// now that the new values of the horizontal velocity values are knoown, the lower boundary condition can be solved
	solve_lower_boundary(state_p2, grid);
	// the advective part of the vertical velocity equation is solved here, the vertical advection is calculated implicitly
	three_band_solver_ver_vel_adv(state_p1, state_p2, state_tendency, delta_t, grid);
	// here, the horizontal divergences are calculated with the new values of the horizontal velocity
	calc_partially_implicit_divvs(state_p2, state_tendency, grid, dualgrid, dissipation_on, 0, tracers_on, delta_t, diffusion_on, radiation_tendency, 0, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divv, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divv, temp_diffusion_heating, temp_gradient, heating_diss, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divv, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divv, 0);
	// here, the non-advective part of the vertical velocity equation is solved implicitly (sound wave solver)
	three_band_solver_ver_sound_waves(state_p1, state_p2, state_tendency, pressure_gradient_acc_1, gradient_geopotential_energy, temperature, temperature_density, temperature_flux_density, temperature_flux_density_divv, wind_field_divv, delta_t, grid);
	// now that the new vertical velocity is known, the new dry density can be calculated via implicit vertical advection
	three_band_solver_ver_den_dry(state_p1, state_p2, state_tendency, delta_t, grid);
	// Here, the global entropy density is calculated via its generalized equation of state.
	three_band_solver_ver_entropy_gas(state_p1, state_p2, state_tendency, delta_t, grid);
	// diagnozing the tendency
	linear_combine_two_states(state_p1, state_p2, state_tendency, -1/delta_t, 1/delta_t);
	// vertical tracer advection with 3-band matrices
	if (tracers_on == 1)
		three_band_solver_ver_tracers(state_p1, state_p2, state_tendency, delta_t, grid);
	// weights of the tendencies
	double first_tendency_weight = 0.0;
	double second_tendency_weight = 1 - first_tendency_weight;
	// dtermining the actual tendency
    linear_combine_two_states(state_tendency_precond, state_tendency, state_tendency, first_tendency_weight, second_tendency_weight);
    // the final step
    linear_combine_two_states(state_0, state_tendency, state_p1, 1, delta_t);
    return 0;
}







