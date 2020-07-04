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

int manage_time_stepping(State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_update_pre, int tracers_on, int diffusion_on, Scalar_field radiation_tendency, double tracer_mass_source_rates[], double tracer_heat_source_rates[], State *state_tendency, Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divv, Scalar_field temperature, Vector_field entropy_gas_flux_density, Scalar_field entropy_gas_flux_density_divv, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Vector_field friction_acc, Scalar_field heating_diss, Scalar_field specific_entropy, Curl_field rel_curl, Curl_field abs_curl, Vector_field downgradient_macroscopic_energy, Vector_field pressure_gradient_acc, Vector_field abs_curl_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field macroscopic_energy, Scalar_field pressure_gradient_decel_factor, Vector_field pressure_gradient_acc_1, Scalar_field diffusion_coeff_numerical_h, Scalar_field diffusion_coeff_numerical_v, Vector_field mass_dry_diffusion_flux_density, Scalar_field mass_dry_diffusion_source_rate, Vector_field temperature_flux_density, Scalar_field tracer_density, Vector_field tracer_velocity, Vector_field tracer_flux_density, Scalar_field tracer_flux_density_divv, Scalar_field tracer_density_temperature, Vector_field tracer_temperature_flux_density, Scalar_field tracer_temperature_flux_density_divv, Vector_field temp_gradient_times_c_h_p)
{
	/*
	Here, the RK3 stepping is implemented.
	*/
	int rad_update, diff_update, diss_update, phase_trans_on;
	for (int i = 0; i < 3; ++i)
		{
		/*
		Diffusion and dissipation are only updated at the first step.
		If radiation is updated, it is done at the first step.
		*/
		if (rad_update_pre == 1 && i == 0)
			rad_update = 1;
		else
			rad_update = 0;
		if (diffusion_on == 1 && i == 0)
			diff_update = 1;
		else
			diff_update = 0;
		if (dissipation_on == 1 && i == 0)
			diss_update = 1;
		else
			diss_update = 0;
		if (tracers_on == 1 && i == 0)
			phase_trans_on = 1;
		else
			phase_trans_on = 0;
		explicit_momentum_tendencies(state_p1, state_tendency, grid, dualgrid, dissipation_on, tracers_on, temperature, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, rel_curl, abs_curl, downgradient_macroscopic_energy, pressure_gradient_acc, abs_curl_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, diss_update, temp_gradient_times_c_h_p);
		three_band_solver_hor(state_0, state_p1, state_tendency, delta_t/(3 - i), grid);
		solve_lower_boundary(state_p1, grid);
		three_band_solver_ver_vel(state_0, state_p1, state_tendency, delta_t/(3 - i), grid);
		calc_partially_implicit_divvs(state_p1, state_tendency, grid, dualgrid, dissipation_on, rad_update, tracers_on, delta_t, diffusion_on, radiation_tendency, phase_trans_on, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divv, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divv, temp_diffusion_heating, temp_gradient, heating_diss, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divv, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divv, diff_update);
		three_band_solver_ver_den_dry(state_0, state_p1, state_tendency, delta_t/(3 - i), grid);
		three_band_solver_ver_entropy_gas(state_0, state_p1, state_tendency, delta_t/(3 - i), grid);
		vertical_sound_wave_solver(state_0, state_p1, state_tendency, delta_t/(3 - i), grid, pressure_gradient_acc, i);
		// implicit vertical tracer advection
		if (i == 2)
			exit(1);
		if (tracers_on == 1)
			three_band_solver_ver_tracers(state_0, state_p1, state_tendency, delta_t/(3 - i), grid);
    }
    return 0;
}












