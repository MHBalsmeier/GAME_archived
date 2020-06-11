#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"
#include "rad.h"
#include "addcomp.h"
#include "surface.h"
#include <stdlib.h>
#include <stdio.h>

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_bool, int add_comps_bool, double delta_t)
{
    Vector_field *density_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density, current_state -> wind, *density_flux, grid);
    Scalar_field *density_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_flux, *density_flux_divergence, grid, 0);
    free(density_flux);
    int retval;
    Vector_field *diffusion_mass_flux = malloc(sizeof(Vector_field));
    Scalar_field *mass_source_rate = malloc(sizeof(Scalar_field));
    Scalar_field *temperature = malloc(sizeof(Scalar_field));
    Scalar_field *mass_diffusion_coeff_numerical_h = malloc(sizeof(Scalar_field));
	Scalar_field *mass_diffusion_coeff_numerical_v = malloc(sizeof(Scalar_field));
	temperature_diagnostics(current_state -> entropy, current_state -> density, current_state -> add_comp_densities, *temperature);
    if (dissipation_on == 1)
    {
        Vector_field *diffusion_mass_flux_pre = malloc(sizeof(Vector_field));
        retval = grad(current_state -> density, *diffusion_mass_flux_pre, grid);
        retval = calc_mass_diffusion_coeffs(*temperature, current_state -> density, *mass_diffusion_coeff_numerical_h, *mass_diffusion_coeff_numerical_v);
        scalar_times_vector_h_v(*mass_diffusion_coeff_numerical_h, *mass_diffusion_coeff_numerical_v, *diffusion_mass_flux_pre, *diffusion_mass_flux, grid);
        free(diffusion_mass_flux_pre);
        retval = divergence(*diffusion_mass_flux, *mass_source_rate, grid, 0);
    }
    free(mass_diffusion_coeff_numerical_h);
	free(mass_diffusion_coeff_numerical_v);
    free(diffusion_mass_flux);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (dissipation_on == 1)
            state_tendency -> density[i] = -(*density_flux_divergence)[i] + (*mass_source_rate)[i];
        else
            state_tendency -> density[i] = -(*density_flux_divergence)[i];
    }
    free(density_flux_divergence);
    Vector_field *entropy_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> entropy, current_state -> wind, *entropy_flux, grid);
    Scalar_field *entropy_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*entropy_flux, *entropy_flux_divergence, grid, 0);
    free(entropy_flux);
    Scalar_field *temp_diffusion_heating = malloc(sizeof(Scalar_field));
    Vector_field *temperature_flux = malloc(sizeof(Vector_field));
    Scalar_field *temp_diffusion_coeff_numerical_h = malloc(sizeof(Scalar_field));
	Scalar_field *temp_diffusion_coeff_numerical_v = malloc(sizeof(Scalar_field));
    if (dissipation_on == 1)
    {
        Vector_field *temperature_flux_pre = malloc(sizeof(Vector_field));
        retval = grad(*temperature, *temperature_flux_pre, grid);
        retval = calc_temp_diffusion_coeffs(*temperature, current_state -> density, *temp_diffusion_coeff_numerical_h, *temp_diffusion_coeff_numerical_v);
        scalar_times_vector_h_v(*temp_diffusion_coeff_numerical_h, *temp_diffusion_coeff_numerical_v, *temperature_flux_pre, *temperature_flux, grid);
        free(temperature_flux_pre);
        retval = divergence(*temperature_flux, *temp_diffusion_heating, grid, 0);
    }
    double rho_h, c_h_v;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	rho_h = current_state -> density[i] + current_state -> add_comp_densities[(NUMBER_OF_ADD_COMPS - 1)*NUMBER_OF_SCALARS + i];
    	c_h_v = spec_heat_cap_diagnostics_v(current_state -> density[i], current_state -> add_comp_densities[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS + i]);
    	(*temp_diffusion_heating)[i] = rho_h*c_h_v*(*temp_diffusion_heating)[i];
    }
    free(temp_diffusion_coeff_numerical_h);
	free(temp_diffusion_coeff_numerical_v);
    free(temperature_flux);
    double viscosity_coeff_molecular = 5*pow(10, -5);
    double viscosity_coeff = pow(10, 5)*viscosity_coeff_molecular;
    Scalar_field *rad_heating = calloc(1, sizeof(Scalar_field));
    if (rad_bool == 1)
        retval = calc_rad_heating(*rad_heating, NUMBER_OF_SCALARS);
    double *add_comp_mass_source_rates = calloc(NUMBER_OF_ADD_COMPS*NUMBER_OF_SCALARS, sizeof(double));
    double *add_comp_heat_source_rates = calloc(NUMBER_OF_ADD_COMPS*NUMBER_OF_SCALARS, sizeof(double));
    if (add_comps_bool == 1)
        retval = calc_add_comp_source_rates(add_comp_mass_source_rates, add_comp_heat_source_rates, current_state -> add_comp_densities, current_state -> add_comp_temps, *temperature, NUMBER_OF_ADD_COMPS, NUMBER_OF_SCALARS, delta_t);
    double total_density;
    Vector_field *friction_acc = malloc(sizeof(Vector_field));
    Scalar_field *heating_diss = malloc(sizeof(Scalar_field));
    retval = dissipation(current_state -> wind, current_state -> density, *friction_acc, *heating_diss, grid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (dissipation_on == 1)
        {
            total_density = current_state -> density[i];
            for (int k = 0; k < NUMBER_OF_ADD_COMPS; ++k)
                total_density += current_state -> add_comp_densities[k*NUMBER_OF_SCALARS + i];
            c_h_v = spec_heat_cap_diagnostics_v(current_state -> density[i], current_state -> add_comp_densities[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS + i]);
            rho_h = current_state -> density[i] + current_state -> add_comp_densities[(NUMBER_OF_ADD_COMPS - 1)*NUMBER_OF_SCALARS + i];
            state_tendency -> entropy[i] = -(*entropy_flux_divergence)[i] + 1/(*temperature)[i]*(rho_h/total_density*((*temp_diffusion_heating)[i] + (*heating_diss)[i] + (*rad_heating)[i]) + add_comp_heat_source_rates[(NUMBER_OF_ADD_COMPS - 1)*NUMBER_OF_SCALARS + i]);
        }
        else
            state_tendency -> entropy[i] = -(*entropy_flux_divergence)[i];
    }
    free(mass_source_rate);
    free(entropy_flux_divergence);
    Scalar_field *add_comp_density = malloc(sizeof(Scalar_field));
    Vector_field *add_comp_velocity = malloc(sizeof(Vector_field));
    Vector_field *add_comp_density_flux = malloc(sizeof(Vector_field));
    Scalar_field *add_comp_flux_divergence = malloc(sizeof(Scalar_field));
    int h_index, layer_index;
    if (add_comps_bool == 1)
    {
		for (int i = 0; i < NUMBER_OF_ADD_COMPS; ++i)
		{
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		        (*add_comp_density)[j] = current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j];
		    if (i < NUMBER_OF_COND_ADD_COMPS)
		    {
		        for (int j = 0; j < NUMBER_OF_VECTORS; ++j)
		        {
		            layer_index = j/NUMBER_OF_VECTORS_PER_LAYER;
		            h_index = j - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
		            (*add_comp_velocity)[j] = current_state -> wind[j];
		            if (h_index < NUMBER_OF_VECTORS_V)
		                (*add_comp_velocity)[j] += ret_sink_velocity(i, 0, 0.001);
		        }
		        retval = scalar_times_vector(*add_comp_density, *add_comp_velocity, *add_comp_density_flux, grid);
		        retval = divergence(*add_comp_density_flux, *add_comp_flux_divergence, grid, 1);
		    }
		    else
		    {
		        retval = scalar_times_vector(*add_comp_density, current_state -> wind, *add_comp_density_flux, grid);
		        retval = divergence(*add_comp_density_flux, *add_comp_flux_divergence, grid, 0);
		    }
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		    {
		        state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] = -(*add_comp_flux_divergence)[j] + add_comp_mass_source_rates[i*NUMBER_OF_SCALARS + j];
		        if (current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j] + delta_t*state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] < 0)
		            state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] = -current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j]/delta_t;
		    }
		}
    }
    else
    {
		for (int i = 0; i < NUMBER_OF_ADD_COMPS; ++i)
		{
		    for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
		        state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] = 0;
		}
    }
    free(add_comp_mass_source_rates);
    free(add_comp_density);
    free(add_comp_density_flux);
    free(add_comp_flux_divergence);
    Scalar_field *add_comp_temp = malloc(sizeof(Scalar_field));
    Scalar_field *add_comp_temp_adv = malloc(sizeof(Scalar_field));
    double c_v_cond;
    for (int i = 0; i < NUMBER_OF_COND_ADD_COMPS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS; ++j)
        {
            layer_index = j/NUMBER_OF_VECTORS_PER_LAYER;
            h_index = j - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
            (*add_comp_velocity)[j] = current_state -> wind[j];
            if (h_index < NUMBER_OF_VECTORS_V)
                (*add_comp_velocity)[j] += ret_sink_velocity(i, 0, 0.001);
        }
        retval = adv_scalar(*add_comp_temp, *add_comp_velocity, *add_comp_temp_adv, grid, dualgrid);
        for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
        {
            c_v_cond = ret_c_v_cond(i, 0, current_state -> add_comp_temps[i*NUMBER_OF_SCALARS + j]);
            total_density = current_state -> density[j];
            for (int k = 0; k < NUMBER_OF_ADD_COMPS; ++k)
                total_density += current_state -> add_comp_densities[k*NUMBER_OF_SCALARS + j];
            if (current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j] > EPSILON_TRACERS)
                state_tendency -> add_comp_temps[i*NUMBER_OF_SCALARS + j] = (*add_comp_temp_adv)[j] + 1/c_v_cond*(1/total_density*((*temp_diffusion_heating)[j] + (*heating_diss)[j] + (*rad_heating)[j]) + 1/current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j]*add_comp_heat_source_rates[i*NUMBER_OF_SCALARS + j]);
            else
                state_tendency -> add_comp_temps[i*NUMBER_OF_SCALARS + j] = -(current_state -> add_comp_temps[i*NUMBER_OF_SCALARS + j] - (*temperature)[j])/delta_t;
        }
    }
    free(temp_diffusion_heating);
    free(heating_diss);
    free(add_comp_heat_source_rates);
    free(rad_heating);
    free(add_comp_temp);
    free(add_comp_velocity);
    free(add_comp_temp_adv);
    Dual_vector_field *rel_curl = malloc(sizeof(Dual_vector_field));
    curl(current_state -> wind, *rel_curl, grid, dualgrid);
    Dual_vector_field *abs_curl = malloc(sizeof(Dual_vector_field));
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        (*abs_curl)[i] = dualgrid -> f_vec[i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER] + (*rel_curl)[i];
    }
    free(rel_curl);
    Scalar_field *specific_entropy = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	rho_h = current_state -> density[i] + current_state -> add_comp_densities[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS + i];
    	(*specific_entropy)[i] = current_state -> entropy[i]/rho_h;
	}
    Vector_field *specific_entropy_gradient = malloc(sizeof(Vector_field));
    retval = grad(*specific_entropy, *specific_entropy_gradient, grid);
    if (retval != 0)
    	printf("grad called at position 0 from tendency errored out, exit code %d.\n", retval);
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
    	(*c_h_p_field)[i] = spec_heat_cap_diagnostics_p(current_state -> density[i], current_state -> add_comp_densities[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS + i]);	
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
    	rho_h = current_state -> density[i] + current_state -> add_comp_densities[(NUMBER_OF_ADD_COMPS - 1)*NUMBER_OF_SCALARS + i];
        total_density = current_state -> density[i];
        for (int k = 0; k < NUMBER_OF_ADD_COMPS; ++k)
            total_density += current_state -> add_comp_densities[k*NUMBER_OF_SCALARS + i];
    	(*pressure_gradient_decel_factor)[i] = rho_h/total_density;
    }
	retval = scalar_times_vector(*pressure_gradient_decel_factor, *pressure_gradient_acc, *pressure_gradient_acc, grid);
	if (retval != 0)
		printf("scalar_times_vector called at position 1 from tendency errored out, exit code %d.\n", retval);
    free(pressure_gradient_decel_factor);
    Vector_field *abs_curl_tend = malloc(sizeof(Vector_field));
    cross_product(current_state -> wind, *abs_curl, *abs_curl_tend, grid);
    free(abs_curl);
    Scalar_field *macroscopic_energy = malloc(sizeof(Scalar_field));
    inner(current_state -> wind, current_state -> wind, *macroscopic_energy, grid, dualgrid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    	(*macroscopic_energy)[i] = -grid -> gravity_potential[i] - 0.5*(*macroscopic_energy)[i];
    Vector_field *downgradient_macroscopic_energy = malloc(sizeof(Vector_field));
    grad(*macroscopic_energy, *downgradient_macroscopic_energy, grid);
    if (retval != 0)
    	printf("grad called at position 2 from tendency errored out, exit code %d.\n", retval);
    free(macroscopic_energy);
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
    	layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (i < NUMBER_OF_VECTORS_V || i >= NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V)
            state_tendency -> wind[i] = 0;
        else
            state_tendency -> wind[i] = (*pressure_gradient_acc)[i] + (*abs_curl_tend)[i] + (*downgradient_macroscopic_energy)[i] + (*friction_acc)[i];
    }
    free(friction_acc);
    free(abs_curl_tend);
    free(pressure_gradient_acc);
    free(downgradient_macroscopic_energy);
    return retval;
}







