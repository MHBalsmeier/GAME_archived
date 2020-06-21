#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int manage_time_stepping(State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_update, int tracers_on, int diffusion_on, Scalar_field radiation_tendency, int tracers_update, int tracers_dynamics_delta_t_ratio, double tracer_mass_source_rates[], double tracer_heat_source_rates[], State *state_tendency, Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divergence, Scalar_field temperature, Vector_field entropy_gas_flux_density, Scalar_field entropy_gas_flux_density_divergence, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Vector_field friction_acc, Scalar_field heating_diss, Scalar_field specific_entropy, Curl_field rel_curl, Curl_field abs_curl, Vector_field downgradient_macroscopic_energy, Vector_field pressure_gradient_acc, Vector_field abs_curl_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field macroscopic_energy, Scalar_field pressure_gradient_decel_factor, Vector_field pressure_gradient_acc_1, Scalar_field diffusion_coeff_numerical_h, Scalar_field diffusion_coeff_numerical_v, Vector_field mass_dry_diffusion_flux_density, Scalar_field mass_dry_diffusion_source_rate, Vector_field temperature_flux_density, Scalar_field tracer_density, Vector_field tracer_velocity, Vector_field tracer_flux_density, Scalar_field tracer_flux_density_divergence, Scalar_field tracer_density_temperature, Vector_field tracer_temperature_flux_density, Scalar_field tracer_temperature_flux_density_divergence)
{
    tendency(state_0, state_tendency, grid, dualgrid, 0, 0, 0, delta_t, 0, radiation_tendency, 0, tracers_dynamics_delta_t_ratio, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divergence, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divergence, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, rel_curl, abs_curl, downgradient_macroscopic_energy, pressure_gradient_acc, abs_curl_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divergence, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divergence);
    State *state_mod_1 = malloc(sizeof(State));
    linear_combine_two_states(state_0, state_tendency, state_mod_1, 1, delta_t);
    State *state_star = malloc(sizeof(State));
    linear_combine_two_states(state_0, state_mod_1, state_star, 2.0/3.0, 1.0/3.0);
    free(state_mod_1);
    tendency(state_star, state_tendency, grid, dualgrid, 0, 0, 0, delta_t, 0, radiation_tendency, 0, tracers_dynamics_delta_t_ratio, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divergence, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divergence, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, rel_curl, abs_curl, downgradient_macroscopic_energy, pressure_gradient_acc, abs_curl_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divergence, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divergence);
    State *state_mod_2 = malloc(sizeof(State));
    linear_combine_two_states(state_0, state_tendency, state_mod_2, 1, delta_t);
    linear_combine_two_states(state_0, state_mod_2, state_star, 0.5, 0.5);
    free(state_mod_2);
    tendency(state_star, state_tendency, grid, dualgrid, dissipation_on, rad_update, tracers_on, delta_t, diffusion_on, radiation_tendency, tracers_on*tracers_update, tracers_dynamics_delta_t_ratio, tracer_mass_source_rates, tracer_heat_source_rates, mass_dry_flux_density, mass_dry_flux_density_divergence, temperature, entropy_gas_flux_density, entropy_gas_flux_density_divergence, temp_diffusion_heating, temp_gradient, friction_acc, heating_diss, specific_entropy, rel_curl, abs_curl, downgradient_macroscopic_energy, pressure_gradient_acc, abs_curl_tend, specific_entropy_gradient, c_h_p_field, macroscopic_energy, pressure_gradient_decel_factor, pressure_gradient_acc_1, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, temperature_flux_density, tracer_density, tracer_velocity, tracer_flux_density, tracer_flux_density_divergence, tracer_density_temperature, tracer_temperature_flux_density, tracer_temperature_flux_density_divergence);
    free(state_star);
    linear_combine_two_states(state_0, state_tendency, state_p1, 1, delta_t);
    int layer_index = NUMBER_OF_LAYERS;
	int layer_index_oro = layer_index - (NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS);
	int retval;
    double u_lowest, v_lowest, n_x, n_y, n_z, check_value;
	for (int h_index = 0; h_index < NUMBER_OF_VECTORS_V; ++h_index)
	{
		n_x = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 0];
		n_y = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 1];
		n_z = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 2];			
		retval = recov_ver_0_pri(state_p1 -> velocity_gas, layer_index, h_index, &u_lowest, grid);
		if (retval != 0)
		{
			printf("Error in recov_ver_pri_0 called at position 0 from manage_time_stepping.\n");
			exit(1);
		}
		retval = recov_ver_1_pri(state_p1 -> velocity_gas, layer_index, h_index, &v_lowest, grid);
		if (retval != 0)
		{
			printf("Error in recov_ver_pri_1 called at position 0 from manage_time_stepping.\n");
			exit(1);
		}
		state_p1 -> velocity_gas[layer_index*NUMBER_OF_VECTORS_PER_LAYER + h_index] = -1/n_z*(n_x*u_lowest + n_y*v_lowest);
		vertical_contravariant_normalized(state_p1 -> velocity_gas, layer_index, h_index, grid, &check_value);
		if (fabs(check_value) > 0.001)
		{
			printf("Error with lower boundary condition.\n");
			exit(1);	
		}
	}
	
    return 0;
}












