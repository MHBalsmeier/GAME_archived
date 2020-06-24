#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include "atmostracers.h"
#include "rte-rrtmgp-c.h"
#include <stdlib.h>
#include <stdio.h>

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_update, int tracers_on, double delta_t, int diffusion_on, Scalar_field radiation_tendency, int tracers_update, int tracers_dynamics_delta_t_ratio, double tracer_mass_source_rates[], double tracer_heat_source_rates[], Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divergence, Scalar_field temperature, Vector_field entropy_gas_flux_density, Scalar_field entropy_gas_flux_density_divergence, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Vector_field friction_acc, Scalar_field heating_diss, Scalar_field specific_entropy, Curl_field rel_curl, Curl_field abs_curl, Vector_field downgradient_macroscopic_energy, Vector_field pressure_gradient_acc, Vector_field abs_curl_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field macroscopic_energy, Scalar_field pressure_gradient_decel_factor, Vector_field pressure_gradient_acc_1, Scalar_field diffusion_coeff_numerical_h, Scalar_field diffusion_coeff_numerical_v, Vector_field mass_dry_diffusion_flux_density, Scalar_field mass_dry_diffusion_source_rate, Vector_field temperature_flux_density, Scalar_field tracer_density, Vector_field tracer_velocity, Vector_field tracer_flux_density, Scalar_field tracer_flux_density_divergence, Scalar_field tracer_density_temperature, Vector_field tracer_temperature_flux_density, Scalar_field tracer_temperature_flux_density_divergence)
{
    scalar_times_vector(current_state -> density_dry, current_state -> velocity_gas, mass_dry_flux_density, grid);
    divergence(mass_dry_flux_density, mass_dry_flux_density_divergence, grid, 0);
	int retval = temperature_diagnostics(current_state -> entropy_gas, current_state -> density_dry, current_state -> tracer_densities, temperature);
    if (diffusion_on == 1)
    {
        retval = grad(current_state -> density_dry, mass_dry_diffusion_flux_density, grid);
		if (retval != 0)
		{
			printf("grad called at position 3 from tendency errored out, exit code %d.\n", retval);
			exit(1);
		}
        retval = calc_mass_diffusion_coeffs(temperature, current_state -> density_dry, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v);
        scalar_times_vector_h_v(diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_flux_density, grid);
        retval = divergence(mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, grid, 0);
		for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	        state_tendency -> density_dry[i] = -mass_dry_flux_density_divergence[i] + mass_dry_diffusion_source_rate[i];
    }
    else
    {
    	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        	state_tendency -> density_dry[i] = -mass_dry_flux_density_divergence[i];
    }
    scalar_times_vector(current_state -> entropy_gas, current_state -> velocity_gas, entropy_gas_flux_density, grid);
    divergence(entropy_gas_flux_density, entropy_gas_flux_density_divergence, grid, 0);
    if (rad_update == 1)
    {
        retval = calc_rad_heating(radiation_tendency, NUMBER_OF_SCALARS);
		if (retval != 0)
		{
			printf("Error in calc_rad_heating called from tendency, position 0.\n");
			exit(1);
		}
    }
    retval = grad(temperature, temp_gradient, grid);
    if (retval != 0)
    {
    	printf("grad called at position 1 from tendency errored out, exit code %d.\n", retval);
    	exit(1);
	}
    double rho_h, c_h_v;
    if (diffusion_on == 1)
    {
        retval = calc_temp_diffusion_coeffs(temperature, current_state -> density_dry, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v);
        scalar_times_vector_h_v(diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, temp_gradient, temperature_flux_density, grid);
        retval = divergence(temperature_flux_density, temp_diffusion_heating, grid, 0);
		for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
		{
			if (tracers_on == 1)
			{
				rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
				c_h_v = spec_heat_cap_diagnostics_v(current_state -> density_dry[i], current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
				temp_diffusion_heating[i] = rho_h*c_h_v*temp_diffusion_heating[i];
			}
			else
				temp_diffusion_heating[i] = current_state -> density_dry[i]*C_D_V*temp_diffusion_heating[i];
		}
    }
    double total_density;
    if (dissipation_on == 1)
    {
		retval = dissipation(current_state -> velocity_gas, current_state -> density_dry, friction_acc, heating_diss, grid);
		if (retval != 0)
		{
			printf("Error in dissipation called from tendency, position 0, exit code was %d.", retval);
			exit(1);
		}
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (diffusion_on == 1)
        {
            total_density = current_state -> density_dry[i];
            for (int k = 0; k < NUMBER_OF_TRACERS; ++k)
                total_density += current_state -> tracer_densities[k*NUMBER_OF_SCALARS + i];
            c_h_v = spec_heat_cap_diagnostics_v(current_state -> density_dry[i], current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
            rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
            state_tendency -> entropy_gas[i] = -entropy_gas_flux_density_divergence[i] + 1/temperature[i]*(rho_h/total_density*(temp_diffusion_heating[i] + heating_diss[i] + radiation_tendency[i]) + tracers_on*tracer_heat_source_rates[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
        }
        else
            state_tendency -> entropy_gas[i] = -entropy_gas_flux_density_divergence[i] + 1/temperature[i]*radiation_tendency[i];
    }
    double c_v_cond;
    int h_index, layer_index;
    if (tracers_update == 1)
    {
        retval = calc_h2otracers_source_rates(tracer_mass_source_rates, tracer_heat_source_rates, current_state -> tracer_densities, current_state -> tracer_density_temperatures, temperature, NUMBER_OF_TRACERS, NUMBER_OF_SCALARS, tracers_dynamics_delta_t_ratio*delta_t);
		if (retval != 0)
		{
			printf("Error in calc_h2otracers_source_rates called from tendency, position 0, exit code was %d.", retval);
			exit(1);
		}
		for (int i = 0; i < NUMBER_OF_TRACERS; ++i)
		{
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		        tracer_density[j] = current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j];
		    if (i < NUMBER_OF_CONDENSATED_TRACERS)
		    {
		        for (int j = 0; j < NUMBER_OF_VECTORS; ++j)
		        {
		            layer_index = j/NUMBER_OF_VECTORS_PER_LAYER;
		            h_index = j - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
		            tracer_velocity[j] = current_state -> velocity_gas[j];
		            if (h_index < NUMBER_OF_VECTORS_V)
		                tracer_velocity[j] -= ret_sink_velocity(i, 0, 0.001);
		        }
		        retval = scalar_times_vector(tracer_density, tracer_velocity, tracer_flux_density, grid);
		        retval = divergence(tracer_flux_density, tracer_flux_density_divergence, grid, 1);
				for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
					tracer_density_temperature[j] = current_state -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j];
		        retval = scalar_times_vector(tracer_density_temperature, tracer_velocity, tracer_temperature_flux_density, grid);
		        retval = divergence(tracer_temperature_flux_density, tracer_temperature_flux_density_divergence, grid, 1);
				for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
				{
				    c_v_cond = ret_c_v_cond(i, 0, current_state -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j]/current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j]);
				    total_density = current_state -> density_dry[j];
				    for (int k = 0; k < NUMBER_OF_TRACERS; ++k)
				        total_density += current_state -> tracer_densities[k*NUMBER_OF_SCALARS + j];
			        state_tendency -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j] = -tracer_temperature_flux_density_divergence[j] + current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j]/(c_v_cond*total_density)*(temp_diffusion_heating[j] + heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*tracer_heat_source_rates[i*NUMBER_OF_SCALARS + j] + tracer_density_temperature[j]*tracer_mass_source_rates[i*NUMBER_OF_SCALARS + j];
			        if (current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j] + tracers_dynamics_delta_t_ratio*delta_t*state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] < 0)
				        state_tendency -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j] = -current_state -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j]/(tracers_dynamics_delta_t_ratio*delta_t);
				}
		    }
		    else
		    {
		        retval = scalar_times_vector(tracer_density, current_state -> velocity_gas, tracer_flux_density, grid);
		        retval = divergence(tracer_flux_density, tracer_flux_density_divergence, grid, 0);
				for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
				{
				    state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] = -tracer_flux_density_divergence[j] + tracer_mass_source_rates[i*NUMBER_OF_SCALARS + j];
				    if (current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j] + tracers_dynamics_delta_t_ratio*delta_t*state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] < 0)
				        state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] = -current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j]/(tracers_dynamics_delta_t_ratio*delta_t);
	            }
		    }
		}
    }
    curl(current_state -> velocity_gas, rel_curl, grid, dualgrid);
    for (int i = 0; i < NUMBER_OF_LAYERS*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H) + NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        layer_index = i/(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
        h_index = i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
        abs_curl[i] = dualgrid -> f_vec[h_index] + rel_curl[i];
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	if (tracers_on == 1)
    	{
			rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
			specific_entropy[i] = current_state -> entropy_gas[i]/rho_h;
    	}
    	else
    		specific_entropy[i] = current_state -> entropy_gas[i]/current_state -> density_dry[i];
	}
    retval = grad(specific_entropy, specific_entropy_gradient, grid);
    if (retval != 0)
    {
    	printf("grad called at position 0 from tendency errored out, exit code %d.\n", retval);
    	exit(1);
	}
	retval = scalar_times_vector(temperature, specific_entropy_gradient, pressure_gradient_acc_1, grid);
	if (retval != 0)
	{
		printf("scalar_times_vector called at position 0 from tendency errored out, exit code %d.\n", retval);
		exit(1);
	}
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	if (tracers_on == 1)
    		c_h_p_field[i] = spec_heat_cap_diagnostics_p(current_state -> density_dry[i], current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
    	else
    		c_h_p_field[i] = C_D_P;
	}
    retval = grad(temperature, temp_gradient, grid);
    if (retval != 0)
    {
    	printf("grad called at position 2 from tendency errored out, exit code %d.\n", retval);
    	exit(1);
	}
	retval = scalar_times_vector(c_h_p_field, temp_gradient, temp_gradient, grid);
	if (retval != 0)
	{
		printf("scalar_times_vector called at position 2 from tendency errored out, exit code %d.\n", retval);
		exit(1);
	}
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
    	pressure_gradient_acc[i] = -temp_gradient[i] + pressure_gradient_acc_1[i];
	}
	if (tracers_on == 1)
	{
		for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
		{
			rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
		    total_density = current_state -> density_dry[i];
		    for (int k = 0; k < NUMBER_OF_TRACERS; ++k)
		    {
		        total_density += current_state -> tracer_densities[k*NUMBER_OF_SCALARS + i];
	        }
			pressure_gradient_decel_factor[i] = rho_h/total_density;
		}
		retval = scalar_times_vector(pressure_gradient_decel_factor, pressure_gradient_acc, pressure_gradient_acc, grid);
		if (retval != 0)
		{
			printf("scalar_times_vector called at position 1 from tendency errored out, exit code %d.\n", retval);
			exit(1);
		}
	}
    retval = coriolis_gen(current_state -> velocity_gas, abs_curl, abs_curl_tend, grid);
    if (retval != 0)
    {
    	printf("Error in coriolis_gen, exit code is %d.\n", retval);
    	exit(1);
    }
    kinetic_energy(current_state -> velocity_gas, macroscopic_energy, grid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    	macroscopic_energy[i] = -grid -> gravity_potential[i] - macroscopic_energy[i];
    grad(macroscopic_energy, downgradient_macroscopic_energy, grid);
    if (retval != 0)
    {
    	printf("grad called at position 4 from tendency errored out, exit code %d.\n", retval);
    	exit(1);
	}
    int layer_index_oro;
    double u_lowest_tendency, v_lowest_tendency, n_x, n_y, n_z;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
    	layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (i < NUMBER_OF_VECTORS_V)
            state_tendency -> velocity_gas[i] = 0;
        else if (i >= NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V)
        {
			layer_index_oro = layer_index - (NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS);
			n_x = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 0];
			n_y = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 1];
			n_z = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 2];			
			retval = recov_ver_0_pri(state_tendency -> velocity_gas, layer_index, h_index, &u_lowest_tendency, grid);
			if (retval != 0)
			{
				printf("Error in recov_ver_pri_0 called at position 0 from horizontal_covariant_normalized.\n");
				exit(1);
			}
			retval = recov_ver_1_pri(state_tendency -> velocity_gas, layer_index, h_index, &v_lowest_tendency, grid);
			if (retval != 0)
			{
				printf("Error in recov_ver_pri_1 called at position 0 from horizontal_covariant_normalized.\n");
				exit(1);
			}
        	state_tendency -> velocity_gas[i] = -1/n_z*(n_x*u_lowest_tendency + n_y*v_lowest_tendency);
        }
        else
        {
        	if (dissipation_on == 1)
            	state_tendency -> velocity_gas[i] = pressure_gradient_acc[i] + abs_curl_tend[i] + downgradient_macroscopic_energy[i] + friction_acc[i];
            else
            	state_tendency -> velocity_gas[i] = pressure_gradient_acc[i] + abs_curl_tend[i] + downgradient_macroscopic_energy[i];
        }
    }
    return retval;
}







