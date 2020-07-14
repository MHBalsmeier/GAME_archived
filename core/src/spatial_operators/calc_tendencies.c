/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include "atmostracers.h"
#include "rte-rrtmgp-c.h"
#include <stdlib.h>
#include <stdio.h>

int explicit_momentum_tendencies(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int tracers_on, Scalar_field temperature, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Vector_field friction_acc, Scalar_field heating_diss, Scalar_field specific_entropy, Curl_field pot_vort, Vector_field gradient_geopotential_energy, Vector_field pressure_gradient_acc, Vector_field pot_vort_tend, Vector_field specific_entropy_gradient, Scalar_field c_h_p_field, Scalar_field e_kin_h, Scalar_field pressure_gradient_decel_factor, Vector_field pressure_gradient_acc_1, int diss_update, Vector_field temp_gradient_times_c_h_p, Vector_field pressure_gradient_acc_old, int no_step, Vector_field e_kin_h_grad, Vector_field mass_dry_flux_density)
{	
	double old_hor_grad_weight = R_D/C_D_P - 0.5;
	double new_hor_grad_weight = C_D_V/C_D_P + 0.5;
	old_hor_grad_weight = 0;
	new_hor_grad_weight = 1;
	temperature_diagnostics(current_state -> entropy_gas, current_state -> density_dry, current_state -> tracer_densities, temperature);
    if (diss_update == 1)
    {
		dissipation(current_state -> velocity_gas, current_state -> density_dry, friction_acc, heating_diss, grid);
    }
    calc_pot_vort(current_state -> velocity_gas, current_state -> density_dry, pot_vort, grid, dualgrid);
    int layer_index, h_index;
    double rho_h;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	if (tracers_on == 1)
    	{
			rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
			specific_entropy[i] = current_state -> entropy_gas[i]/rho_h;
    	}
    	else
    		specific_entropy[i] = current_state -> entropy_gas[i]/current_state -> density_dry[i];
	}
    grad(specific_entropy, specific_entropy_gradient, grid);
	scalar_times_vector(temperature, specific_entropy_gradient, pressure_gradient_acc_1, grid);
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	if (tracers_on == 1)
    		c_h_p_field[i] = spec_heat_cap_diagnostics_p(current_state -> density_dry[i], current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
    	else
    		c_h_p_field[i] = C_D_P;
	}
    grad(temperature, temp_gradient, grid);
	scalar_times_vector(c_h_p_field, temp_gradient, temp_gradient_times_c_h_p, grid);
	if (no_step == 0)
	{
		*pressure_gradient_acc_old = *pressure_gradient_acc;
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			pressure_gradient_acc[i] = -temp_gradient_times_c_h_p[i] + pressure_gradient_acc_1[i];
		}
	}
	else
	{
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			layer_index = i/NO_OF_VECTORS_PER_LAYER;
			h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
			if (h_index < NO_OF_VECTORS_V)
				pressure_gradient_acc[i] = -temp_gradient_times_c_h_p[i] + pressure_gradient_acc_1[i];
		}	
	}
	double total_density;
	if (tracers_on == 1)
	{
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
		    total_density = current_state -> density_dry[i];
		    for (int k = 0; k < NO_OF_TRACERS; ++k)
		    {
		        total_density += current_state -> tracer_densities[k*NO_OF_SCALARS + i];
	        }
			pressure_gradient_decel_factor[i] = rho_h/total_density;
		}
		scalar_times_vector(pressure_gradient_decel_factor, pressure_gradient_acc, pressure_gradient_acc, grid);
	}
    scalar_times_vector(current_state -> density_dry, current_state -> velocity_gas, mass_dry_flux_density, grid);
    coriolis_gen(mass_dry_flux_density, pot_vort, pot_vort_tend, grid);
    kinetic_energy(current_state -> velocity_gas, e_kin_h, grid);
    grad(e_kin_h, e_kin_h_grad, grid);
    grad(grid -> gravity_potential, gradient_geopotential_energy, grid);
    int hor_non_trad_cori_sign;
    double metric_term, vertical_velocity, hor_non_trad_cori_term;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (i < NO_OF_VECTORS_V || i >= NO_OF_VECTORS - NO_OF_VECTORS_V)
            state_tendency -> velocity_gas[i] = 0;
        else
        {
        	if (dissipation_on == 1)
        	{
        		if (h_index >= NO_OF_VECTORS_V)
        		{
        			recov_hor_ver_pri(current_state -> velocity_gas, layer_index, h_index - NO_OF_VECTORS_V, &vertical_velocity, grid);
        			metric_term = -vertical_velocity*current_state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
        			hor_non_trad_cori_sign = -dualgrid -> h_curl_signs[4*(h_index - NO_OF_VECTORS_V) + 0];
        			hor_non_trad_cori_term = hor_non_trad_cori_sign*vertical_velocity*dualgrid -> f_vec[h_index - NO_OF_VECTORS_V];
            		state_tendency -> velocity_gas[i] = old_hor_grad_weight*pressure_gradient_acc_old[i] + new_hor_grad_weight*pressure_gradient_acc[i] + pot_vort_tend[i] - gradient_geopotential_energy[i] - e_kin_h_grad[i] + hor_non_trad_cori_term + metric_term + friction_acc[i];
        		}
            	if (h_index < NO_OF_VECTORS_V)
            	{
            		state_tendency -> velocity_gas[i] = pot_vort_tend[i] - e_kin_h_grad[i] + friction_acc[i];
        		}
        	}
            else
            {
        		if (h_index >= NO_OF_VECTORS_V)
        		{
        			recov_hor_ver_pri(current_state -> velocity_gas, layer_index, h_index - NO_OF_VECTORS_V, &vertical_velocity, grid);
        			metric_term = -vertical_velocity*current_state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
        			hor_non_trad_cori_sign = -dualgrid -> h_curl_signs[4*(h_index - NO_OF_VECTORS_V) + 0];
        			hor_non_trad_cori_term = hor_non_trad_cori_sign*vertical_velocity*dualgrid -> f_vec[h_index - NO_OF_VECTORS_V];
            		state_tendency -> velocity_gas[i] = old_hor_grad_weight*pressure_gradient_acc_old[i] + new_hor_grad_weight*pressure_gradient_acc[i] + pot_vort_tend[i] - gradient_geopotential_energy[i] - e_kin_h_grad[i] + hor_non_trad_cori_term + metric_term;
        		}
        		if (h_index < NO_OF_VECTORS_V)
        		{
            		state_tendency -> velocity_gas[i] = pot_vort_tend[i] - e_kin_h_grad[i];
        		}
        	}
        }
    }
    return 0;
}

int calc_partially_implicit_divvs(State *next_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_update, int tracers_on, double delta_t, int diffusion_on, Scalar_field radiation_tendency, int phase_transitions_on, double tracer_mass_source_rates[], double tracer_heat_source_rates[], Vector_field mass_dry_flux_density, Scalar_field mass_dry_flux_density_divv, Scalar_field temperature, Vector_field entropy_gas_flux_density, Scalar_field entropy_gas_flux_density_divv, Scalar_field temp_diffusion_heating, Vector_field temp_gradient, Scalar_field heating_diss, Scalar_field diffusion_coeff_numerical_h, Scalar_field diffusion_coeff_numerical_v, Vector_field mass_dry_diffusion_flux_density, Scalar_field mass_dry_diffusion_source_rate, Vector_field temperature_flux_density, Scalar_field tracer_density, Vector_field tracer_velocity, Vector_field tracer_flux_density, Scalar_field tracer_flux_density_divv, Scalar_field tracer_density_temperature, Vector_field tracer_temperature_flux_density, Scalar_field tracer_temperature_flux_density_divv, int diff_update)
{
    scalar_times_vector(next_state -> density_dry, next_state -> velocity_gas, mass_dry_flux_density, grid);
    divv_h(mass_dry_flux_density, mass_dry_flux_density_divv, grid);
    if (diff_update == 1)
    {
        grad(next_state -> density_dry, mass_dry_diffusion_flux_density, grid);
        calc_mass_diffusion_coeffs(temperature, next_state -> density_dry, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v);
        scalar_times_vector_scalar_h_v(diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, mass_dry_diffusion_flux_density, mass_dry_diffusion_flux_density, grid);
        divv(mass_dry_diffusion_flux_density, mass_dry_diffusion_source_rate, grid, 0);
    }
    if (diffusion_on == 1)
    {
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    state_tendency -> density_dry[i] = -mass_dry_flux_density_divv[i] + mass_dry_diffusion_source_rate[i];
	    }
    }
    else
    {
    	for (int i = 0; i < NO_OF_SCALARS; ++i)
        	state_tendency -> density_dry[i] = -mass_dry_flux_density_divv[i];
    }
    scalar_times_vector(next_state -> entropy_gas, next_state -> velocity_gas, entropy_gas_flux_density, grid);
    divv_h(entropy_gas_flux_density, entropy_gas_flux_density_divv, grid);
    if (rad_update == 1)
    {
        calc_rad_heating(radiation_tendency, NO_OF_SCALARS);
    }
    double rho_h, c_h_v;
    if (diff_update == 1)
    {
        calc_temp_diffusion_coeffs(temperature, next_state -> density_dry, diffusion_coeff_numerical_h, diffusion_coeff_numerical_v);
        scalar_times_vector_scalar_h_v(diffusion_coeff_numerical_h, diffusion_coeff_numerical_v, temp_gradient, temperature_flux_density, grid);
        divv(temperature_flux_density, temp_diffusion_heating, grid, 0);
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (tracers_on == 1)
			{
				rho_h = next_state -> density_dry[i] + next_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
				c_h_v = spec_heat_cap_diagnostics_v(next_state -> density_dry[i], next_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
				temp_diffusion_heating[i] = rho_h*c_h_v*temp_diffusion_heating[i];
			}
			else
				temp_diffusion_heating[i] = next_state -> density_dry[i]*C_D_V*temp_diffusion_heating[i];
		}
    }
    double total_density;
    double c_v_cond;
    int h_index, layer_index;
    if (tracers_on == 1)
    {
    	/*
    	phase transitions are on only at the third RK step
    	only then, they are also updated
    	*/
    	if (phase_transitions_on == 1)
    	{
		    calc_h2otracers_source_rates(tracer_mass_source_rates, tracer_heat_source_rates, next_state -> tracer_densities, next_state -> tracer_density_temperatures, temperature, NO_OF_TRACERS, NO_OF_SCALARS, delta_t);
		}
		for (int i = 0; i < NO_OF_TRACERS; ++i)
		{
		    for (int j = 0; j < NO_OF_SCALARS; ++j)
		        tracer_density[j] = next_state -> tracer_densities[i*NO_OF_SCALARS + j];
		    if (i < NO_OF_CONDENSATED_TRACERS)
		    {
		        for (int j = 0; j < NO_OF_VECTORS; ++j)
		        {
		            layer_index = j/NO_OF_VECTORS_PER_LAYER;
		            h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
		            if (h_index < NO_OF_VECTORS_V)
		                tracer_velocity[j] = next_state -> velocity_gas[j] - ret_sink_velocity(i, 0, 0.001);
		            else
		            	tracer_velocity[j] = next_state -> velocity_gas[j];
		        }
		        scalar_times_vector(tracer_density, tracer_velocity, tracer_flux_density, grid);
		        divv_h(tracer_flux_density, tracer_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
					tracer_density_temperature[j] = next_state -> tracer_density_temperatures[i*NO_OF_SCALARS + j];
		        scalar_times_vector(tracer_density_temperature, tracer_velocity, tracer_temperature_flux_density, grid);
		        divv_h(tracer_temperature_flux_density, tracer_temperature_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
				    c_v_cond = ret_c_v_cond(i, 0, next_state -> tracer_density_temperatures[i*NO_OF_SCALARS + j]/next_state -> tracer_densities[i*NO_OF_SCALARS + j]);
				    total_density = next_state -> density_dry[j];
				    for (int k = 0; k < NO_OF_TRACERS; ++k)
				        total_density += next_state -> tracer_densities[k*NO_OF_SCALARS + j];
			        state_tendency -> tracer_density_temperatures[i*NO_OF_SCALARS + j] = -tracer_temperature_flux_density_divv[j] + next_state -> tracer_densities[i*NO_OF_SCALARS + j]/(c_v_cond*total_density)*(temp_diffusion_heating[j] + heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*phase_transitions_on*tracer_heat_source_rates[i*NO_OF_SCALARS + j] + tracer_density_temperature[j]*phase_transitions_on*tracer_mass_source_rates[i*NO_OF_SCALARS + j];
				}
		    }
		    else
		    {
		        scalar_times_vector_vector_h_v(tracer_density, next_state -> velocity_gas, next_state -> velocity_gas, tracer_flux_density, grid);
		        divv_h(tracer_flux_density, tracer_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
				    state_tendency -> tracer_densities[i*NO_OF_SCALARS + j] = -tracer_flux_density_divv[j] + phase_transitions_on*tracer_mass_source_rates[i*NO_OF_SCALARS + j];
	            }
		    }
		}
    }
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        if (diffusion_on == 1)
        {
            total_density = next_state -> density_dry[i];
            for (int k = 0; k < NO_OF_TRACERS; ++k)
                total_density += next_state -> tracer_densities[k*NO_OF_SCALARS + i];
            c_h_v = spec_heat_cap_diagnostics_v(next_state -> density_dry[i], next_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
            rho_h = next_state -> density_dry[i] + next_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
            state_tendency -> entropy_gas[i] = -entropy_gas_flux_density_divv[i] + 1/temperature[i]*(rho_h/total_density*(temp_diffusion_heating[i] + dissipation_on*heating_diss[i] + radiation_tendency[i]) + tracers_on*phase_transitions_on*tracer_heat_source_rates[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
        }
        else
        {
            state_tendency -> entropy_gas[i] = -entropy_gas_flux_density_divv[i] + 1/temperature[i]*radiation_tendency[i];
        }
    }
    return 0;
}


















