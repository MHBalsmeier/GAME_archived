#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include "atmosrad.h"
#include "atmostracers.h"
#include "surface.h"
#include <stdlib.h>
#include <stdio.h>

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_on, int tracers_on, double delta_t, int diffusion_on)
{
    Vector_field *mass_dry_flux_density = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density_dry, current_state -> velocity_gas, *mass_dry_flux_density, grid);
    Scalar_field *mass_dry_flux_density_divergence = malloc(sizeof(Scalar_field));
    divergence(*mass_dry_flux_density, *mass_dry_flux_density_divergence, grid, 0);
    free(mass_dry_flux_density);
    Vector_field *mass_dry_diffusion_flux_density = malloc(sizeof(Vector_field));
    Scalar_field *mass_dry_diffusion_source_rate = malloc(sizeof(Scalar_field));
    Scalar_field *temperature = malloc(sizeof(Scalar_field));
    Scalar_field *mass_diffusion_coeff_numerical_h = malloc(sizeof(Scalar_field));
	Scalar_field *mass_diffusion_coeff_numerical_v = malloc(sizeof(Scalar_field));
	int retval = temperature_diagnostics(current_state -> entropy_gas, current_state -> density_dry, current_state -> tracer_densities, *temperature);
    if (diffusion_on == 1)
    {
        Vector_field *mass_dry_diffusion_flux_density_pre = malloc(sizeof(Vector_field));
        retval = grad(current_state -> density_dry, *mass_dry_diffusion_flux_density_pre, grid);
        retval = calc_mass_diffusion_coeffs(*temperature, current_state -> density_dry, *mass_diffusion_coeff_numerical_h, *mass_diffusion_coeff_numerical_v);
        scalar_times_vector_h_v(*mass_diffusion_coeff_numerical_h, *mass_diffusion_coeff_numerical_v, *mass_dry_diffusion_flux_density_pre, *mass_dry_diffusion_flux_density, grid);
        free(mass_dry_diffusion_flux_density_pre);
        retval = divergence(*mass_dry_diffusion_flux_density, *mass_dry_diffusion_source_rate, grid, 0);
    }
    free(mass_diffusion_coeff_numerical_h);
	free(mass_diffusion_coeff_numerical_v);
    free(mass_dry_diffusion_flux_density);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (diffusion_on == 1)
            state_tendency -> density_dry[i] = -(*mass_dry_flux_density_divergence)[i] + (*mass_dry_diffusion_source_rate)[i];
        else
            state_tendency -> density_dry[i] = -(*mass_dry_flux_density_divergence)[i];
    }
    free(mass_dry_flux_density_divergence);
    free(mass_dry_diffusion_source_rate);
    Vector_field *entropy_gas_flux_density = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> entropy_gas, current_state -> velocity_gas, *entropy_gas_flux_density, grid);
    Scalar_field *entropy_gas_flux_density_divergence = malloc(sizeof(Scalar_field));
    divergence(*entropy_gas_flux_density, *entropy_gas_flux_density_divergence, grid, 0);
    free(entropy_gas_flux_density);
    Scalar_field *temp_diffusion_heating = malloc(sizeof(Scalar_field));
    Vector_field *temperature_flux_density = malloc(sizeof(Vector_field));
    Scalar_field *temp_diffusion_coeff_numerical_h = malloc(sizeof(Scalar_field));
	Scalar_field *temp_diffusion_coeff_numerical_v = malloc(sizeof(Scalar_field));
    if (diffusion_on == 1)
    {
        Vector_field *temperature_flux_density_pre = malloc(sizeof(Vector_field));
        retval = grad(*temperature, *temperature_flux_density_pre, grid);
        retval = calc_temp_diffusion_coeffs(*temperature, current_state -> density_dry, *temp_diffusion_coeff_numerical_h, *temp_diffusion_coeff_numerical_v);
        scalar_times_vector_h_v(*temp_diffusion_coeff_numerical_h, *temp_diffusion_coeff_numerical_v, *temperature_flux_density_pre, *temperature_flux_density, grid);
        free(temperature_flux_density_pre);
        retval = divergence(*temperature_flux_density, *temp_diffusion_heating, grid, 0);
    }
    double rho_h, c_h_v;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
    	c_h_v = spec_heat_cap_diagnostics_v(current_state -> density_dry[i], current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
    	(*temp_diffusion_heating)[i] = rho_h*c_h_v*(*temp_diffusion_heating)[i];
    }
    free(temp_diffusion_coeff_numerical_h);
	free(temp_diffusion_coeff_numerical_v);
    free(temperature_flux_density);
    Scalar_field *rad_heating = calloc(1, sizeof(Scalar_field));
    if (rad_on == 1)
        retval = calc_rad_heating(*rad_heating, NUMBER_OF_SCALARS);
    double *tracer_mass_source_rates = calloc(NUMBER_OF_TRACERS*NUMBER_OF_SCALARS, sizeof(double));
    double *tracer_heat_source_rates = calloc(NUMBER_OF_TRACERS*NUMBER_OF_SCALARS, sizeof(double));
    if (tracers_on == 1)
    {
        retval = calc_h2otracers_source_rates(tracer_mass_source_rates, tracer_heat_source_rates, current_state -> tracer_densities, current_state -> tracer_density_temperatures, *temperature, NUMBER_OF_TRACERS, NUMBER_OF_SCALARS, delta_t);
		if (retval != 0)
		{
			printf("Error in calc_h2otracers_source_rates called from tendency, position 0, exit code was %d.", retval);
			exit(1);
		}
    }
    double total_density;
    Vector_field *friction_acc = malloc(sizeof(Vector_field));
    Scalar_field *heating_diss = malloc(sizeof(Scalar_field));
    if (dissipation_on == 1)
    {
		retval = dissipation(current_state -> velocity_gas, current_state -> density_dry, *friction_acc, *heating_diss, grid);
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
            state_tendency -> entropy_gas[i] = -(*entropy_gas_flux_density_divergence)[i] + 1/(*temperature)[i]*(rho_h/total_density*((*temp_diffusion_heating)[i] + (*heating_diss)[i] + (*rad_heating)[i]) + tracers_on*tracer_heat_source_rates[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
        }
        else
            state_tendency -> entropy_gas[i] = -(*entropy_gas_flux_density_divergence)[i];
    }
    free(entropy_gas_flux_density_divergence);
    double c_v_cond;
    int h_index, layer_index;
    if (tracers_on == 1)
    {
		Scalar_field *tracer_density = malloc(sizeof(Scalar_field));
		Vector_field *tracer_velocity = malloc(sizeof(Vector_field));
		Vector_field *tracer_flux_density = malloc(sizeof(Vector_field));
		Scalar_field *tracer_flux_density_divergence = malloc(sizeof(Scalar_field));
		Scalar_field *tracer_density_temperature = malloc(sizeof(Scalar_field));
		Vector_field *tracer_temperature_flux_density = malloc(sizeof(Vector_field));
		Scalar_field *tracer_temperature_flux_density_divergence = malloc(sizeof(Scalar_field));
		for (int i = 0; i < NUMBER_OF_TRACERS; ++i)
		{
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		        (*tracer_density)[j] = current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j];
		    if (i < NUMBER_OF_CONDENSATED_TRACERS)
		    {
		        for (int j = 0; j < NUMBER_OF_VECTORS; ++j)
		        {
		            layer_index = j/NUMBER_OF_VECTORS_PER_LAYER;
		            h_index = j - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
		            (*tracer_velocity)[j] = current_state -> velocity_gas[j];
		            if (h_index < NUMBER_OF_VECTORS_V)
		                (*tracer_velocity)[j] -= ret_sink_velocity(i, 0, 0.001);
		        }
		        retval = scalar_times_vector(*tracer_density, *tracer_velocity, *tracer_flux_density, grid);
		        retval = divergence(*tracer_flux_density, *tracer_flux_density_divergence, grid, 1);
				for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
					(*tracer_density_temperature)[j] = current_state -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j];
		        retval = scalar_times_vector(*tracer_density_temperature, *tracer_velocity, *tracer_temperature_flux_density, grid);
		        retval = divergence(*tracer_temperature_flux_density, *tracer_temperature_flux_density_divergence, grid, 1);
				for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
				{
				    c_v_cond = ret_c_v_cond(i, 0, current_state -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j]/current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j]);
				    total_density = current_state -> density_dry[j];
				    for (int k = 0; k < NUMBER_OF_TRACERS; ++k)
				        total_density += current_state -> tracer_densities[k*NUMBER_OF_SCALARS + j];
			        state_tendency -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j] = -(*tracer_temperature_flux_density_divergence)[j] + current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j]/(c_v_cond*total_density)*((*temp_diffusion_heating)[j] + (*heating_diss)[j] + (*rad_heating)[j]) + 1/c_v_cond*tracer_heat_source_rates[i*NUMBER_OF_SCALARS + j] + (*tracer_density_temperature)[j]*tracer_mass_source_rates[i*NUMBER_OF_SCALARS + j];
			        if (current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j] + delta_t*state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] < 0)
				        state_tendency -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j] = -current_state -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j]/delta_t;
				}
		    }
		    else
		    {
		        retval = scalar_times_vector(*tracer_density, current_state -> velocity_gas, *tracer_flux_density, grid);
		        retval = divergence(*tracer_flux_density, *tracer_flux_density_divergence, grid, 0);
				for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
				{
				    state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] = -(*tracer_flux_density_divergence)[j] + tracer_mass_source_rates[i*NUMBER_OF_SCALARS + j];
				    if (current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j] + delta_t*state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] < 0)
				        state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] = -current_state -> tracer_densities[i*NUMBER_OF_SCALARS + j]/delta_t;
	            }
		    }
		}
		free(tracer_density);
		free(tracer_velocity);
		free(tracer_flux_density);
		free(tracer_flux_density_divergence);
		free(tracer_density_temperature);
		free(tracer_temperature_flux_density);
		free(tracer_temperature_flux_density_divergence);
    }
    else
    {
		for (int i = 0; i < NUMBER_OF_TRACERS; ++i)
		{
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		        state_tendency -> tracer_densities[i*NUMBER_OF_SCALARS + j] = 0;
		}
		for (int i = 0; i < NUMBER_OF_CONDENSATED_TRACERS; ++i)
		{
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		        state_tendency -> tracer_density_temperatures[i*NUMBER_OF_SCALARS + j] = 0;
		}
    }
    free(tracer_mass_source_rates);
    free(temp_diffusion_heating);
    free(tracer_heat_source_rates);
    free(heating_diss);
    free(rad_heating);
    Curl_field *rel_curl = malloc(sizeof(Curl_field));
    curl(current_state -> velocity_gas, *rel_curl, grid, dualgrid);
    Curl_field *abs_curl = malloc(sizeof(Curl_field));
    for (int i = 0; i < NUMBER_OF_LAYERS*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H) + NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        layer_index = i/(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
        h_index = i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
        (*abs_curl)[i] = dualgrid -> f_vec[h_index] + (*rel_curl)[i];
    }
    free(rel_curl);
    Scalar_field *specific_entropy = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
    	(*specific_entropy)[i] = current_state -> entropy_gas[i]/rho_h;
	}
    Vector_field *specific_entropy_gradient = malloc(sizeof(Vector_field));
    retval = grad(*specific_entropy, *specific_entropy_gradient, grid);
    if (retval != 0)
    {
    	printf("grad called at position 0 from tendency errored out, exit code %d.\n", retval);
	}
	free(specific_entropy);
	Vector_field *pressure_gradient_acc_1 = malloc(sizeof(Vector_field));
	retval = scalar_times_vector(*temperature, *specific_entropy_gradient, *pressure_gradient_acc_1, grid);
	if (retval != 0)
		printf("scalar_times_vector called at position 0 from tendency errored out, exit code %d.\n", retval);
	free(specific_entropy_gradient);
	Vector_field *temp_gradient_times_c_h_p = malloc(sizeof(Vector_field));
    retval = grad(*temperature, *temp_gradient_times_c_h_p, grid);
    if (retval != 0)
    	printf("grad called at position 1 from tendency errored out, exit code %d.\n", retval);
    free(temperature);
    Scalar_field *c_h_p_field = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	(*c_h_p_field)[i] = spec_heat_cap_diagnostics_p(current_state -> density_dry[i], current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i]);
	}
	retval = scalar_times_vector(*c_h_p_field, *temp_gradient_times_c_h_p, *temp_gradient_times_c_h_p, grid);
	if (retval != 0)
		printf("scalar_times_vector called at position 2 from tendency errored out, exit code %d.\n", retval);
	free(c_h_p_field);
    Vector_field *pressure_gradient_acc = malloc(sizeof(Vector_field));
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    	(*pressure_gradient_acc)[i] = -(*temp_gradient_times_c_h_p)[i] + (*pressure_gradient_acc_1)[i];
    free(pressure_gradient_acc_1);
    free(temp_gradient_times_c_h_p);
    Scalar_field *pressure_gradient_decel_factor = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i];
        total_density = current_state -> density_dry[i];
        for (int k = 0; k < NUMBER_OF_TRACERS; ++k)
            total_density += current_state -> tracer_densities[k*NUMBER_OF_SCALARS + i];
    	(*pressure_gradient_decel_factor)[i] = rho_h/total_density;
    }
	retval = scalar_times_vector(*pressure_gradient_decel_factor, *pressure_gradient_acc, *pressure_gradient_acc, grid);
	if (retval != 0)
	{
		printf("scalar_times_vector called at position 1 from tendency errored out, exit code %d.\n", retval);
		exit(1);
	}
    free(pressure_gradient_decel_factor);
    Vector_field *abs_curl_tend = malloc(sizeof(Vector_field));
    retval = coriolis_gen(current_state -> velocity_gas, *abs_curl, *abs_curl_tend, grid);
    if (retval != 0)
    {
    	printf(".Error in coriolis_gen, exit code is %d.\n", retval);
    	exit(1);
    }
    free(abs_curl);
    Scalar_field *macroscopic_energy = malloc(sizeof(Scalar_field));
    kinetic_energy(current_state -> velocity_gas, *macroscopic_energy, grid, dualgrid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    	(*macroscopic_energy)[i] = -grid -> gravity_potential[i] - (*macroscopic_energy)[i];
    Vector_field *downgradient_macroscopic_energy = malloc(sizeof(Vector_field));
    grad(*macroscopic_energy, *downgradient_macroscopic_energy, grid);
    if (retval != 0)
    	printf("grad called at position 2 from tendency errored out, exit code %d.\n", retval);
    free(macroscopic_energy);
    int layer_index_oro;
    double u_lowest_tendency, v_lowest_tendency, n_x, n_y, n_z;
    double check_value;
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
				printf("Error in recov_ver_pri_0 called at position 0 from horizontal_covariant_normalized.\n");
			retval = recov_ver_1_pri(state_tendency -> velocity_gas, layer_index, h_index, &v_lowest_tendency, grid);
			if (retval != 0)
				printf("Error in recov_ver_pri_1 called at position 0 from horizontal_covariant_normalized.\n");
        	state_tendency -> velocity_gas[i] = -1/n_z*(n_x*u_lowest_tendency + n_y*v_lowest_tendency);
        	vertical_contravariant_normalized(current_state -> velocity_gas, layer_index, h_index, grid, &check_value);
        	if (fabs(check_value) > 0.001)
        	{
        		printf("Error with lower boundary condition.\n");
        		exit(1);	
    		}
        }
        else
        {
        	if (dissipation_on == 1)
            	state_tendency -> velocity_gas[i] = (*pressure_gradient_acc)[i] + (*abs_curl_tend)[i] + (*downgradient_macroscopic_energy)[i] + (*friction_acc)[i];
            else
            	state_tendency -> velocity_gas[i] = (*pressure_gradient_acc)[i] + (*abs_curl_tend)[i] + (*downgradient_macroscopic_energy)[i];
        }
    }
    free(friction_acc);
    free(abs_curl_tend);
    free(pressure_gradient_acc);
    free(downgradient_macroscopic_energy);
    return retval;
}







