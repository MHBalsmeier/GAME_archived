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

int explicit_momentum_tendencies(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	// Here, weights of the horizontal pressure gradient can be chosen as usual.
	// This is in potential temperature formulation.
	// RES_ID = 5
	// R_D/C_D_P - 0.5:  hrs to crash
	// -2.0/3: 12.09 hrs to crash
	// -0.7:	12.35 hrs to crash
	// -0.8:	13.09 hrs to crash
	// -0.9:	13.09 hrs to crash
	// -1.0:	13.09 hrs to crash
	double old_hor_grad_weight = -0.5;
	double new_hor_grad_weight = 1 - old_hor_grad_weight;
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
	}
	// old_hor_grad_weight = 0;
	// new_hor_grad_weight = 1;
	temperature_diagnostics(current_state, diagnostics -> temperature);
    if (config_info -> momentum_diffusion_on == 1)
    {
		dissipation(current_state -> velocity_gas, current_state -> density_dry, diffusion_info -> friction_acc, diffusion_info -> heating_diss, grid);
    }
    calc_pot_vort(current_state -> velocity_gas, current_state -> density_dry, diagnostics -> pot_vort, grid, dualgrid);
    int layer_index, h_index;
    double rho_h;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	if (config_info -> tracers_on == 1)
    	{
			rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
			diagnostics -> specific_entropy[i] = current_state -> entropy_gas[i]/rho_h;
    	}
    	else
    		diagnostics -> specific_entropy[i] = current_state -> entropy_gas[i]/current_state -> density_dry[i];
	}
    grad(diagnostics -> specific_entropy, diagnostics -> specific_entropy_gradient, grid);
	scalar_times_vector(diagnostics -> temperature, diagnostics -> specific_entropy_gradient, diagnostics -> pressure_gradient_acc_1, grid);
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	if (config_info -> tracers_on == 1)
    		diagnostics -> c_h_p_field[i] = spec_heat_cap_diagnostics_p(current_state -> density_dry[i], current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
    	else
    		diagnostics -> c_h_p_field[i] = C_D_P;
	}
    grad(diagnostics -> temperature, diagnostics -> temp_gradient, grid);
	scalar_times_vector(diagnostics -> c_h_p_field, diagnostics -> temp_gradient, forcings -> temp_gradient_times_c_h_p, grid);
	// Here, the update of the pressure gradient is managed.
	if (no_step_rk == 0)
	{
		if (config_info -> totally_first_step_bool == 0)
		{
			for (int i = 0; i < NO_OF_VECTORS; ++i)
			{
				interpolation -> pressure_gradient_acc_old[i] = forcings -> pressure_gradient_acc[i];
			}
		}
		else
		{
			for (int i = 0; i < NO_OF_VECTORS; ++i)
			{
				interpolation -> pressure_gradient_acc_old[i] = 0;
			}
		}
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			forcings -> pressure_gradient_acc[i] = -forcings -> temp_gradient_times_c_h_p[i] + diagnostics -> pressure_gradient_acc_1[i];
		}
	}
	else
	{
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			layer_index = i/NO_OF_VECTORS_PER_LAYER;
			h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
			if (h_index < NO_OF_VECTORS_V)
				forcings -> pressure_gradient_acc[i] = -forcings -> temp_gradient_times_c_h_p[i] + diagnostics -> pressure_gradient_acc_1[i];
		}	
	}
	// The pressure gradient has to get a deceleration factor in presence of condensates.
	double total_density;
	if (config_info -> tracers_on == 1)
	{
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			rho_h = current_state -> density_dry[i] + current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
		    total_density = current_state -> density_dry[i];
		    for (int k = 0; k < NO_OF_TRACERS; ++k)
		    {
		        total_density += current_state -> tracer_densities[k*NO_OF_SCALARS + i];
	        }
			diffusion_info -> pressure_gradient_decel_factor[i] = rho_h/total_density;
		}
		scalar_times_vector(diffusion_info -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc, forcings -> pressure_gradient_acc, grid);
	}
	// Here, the gaseous flux density is prepared fore the generalized Coriolis term.
    scalar_times_vector(current_state -> density_dry, current_state -> velocity_gas, diagnostics -> mass_dry_flux_density, grid);
    // Noew, the generalized Coriolis term is evaluated.
    coriolis_gen(diagnostics -> mass_dry_flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid);
    // Horizontal kinetic energy is prepared for the gradient term of the Lamb transformation.
    kinetic_energy(current_state -> velocity_gas, diagnostics -> e_kin_h, grid);
    grad(diagnostics -> e_kin_h, forcings -> e_kin_h_grad, grid);
    // Now the explicit forces are added up.
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
        	if (config_info -> momentum_diffusion_on == 1)
        	{
        		if (h_index >= NO_OF_VECTORS_V)
        		{
        			recov_hor_ver_pri(current_state -> velocity_gas, layer_index, h_index - NO_OF_VECTORS_V, &vertical_velocity, grid);
        			metric_term = -vertical_velocity*current_state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
        			hor_non_trad_cori_sign = -dualgrid -> h_curl_signs[4*(h_index - NO_OF_VECTORS_V) + 0];
        			hor_non_trad_cori_term = hor_non_trad_cori_sign*vertical_velocity*dualgrid -> f_vec[h_index - NO_OF_VECTORS_V];
            		state_tendency -> velocity_gas[i] = old_hor_grad_weight*interpolation -> pressure_gradient_acc_old[i] + new_hor_grad_weight*forcings -> pressure_gradient_acc[i] + forcings -> pot_vort_tend[i] - grid -> gravity_m[i] - forcings -> e_kin_h_grad[i] + hor_non_trad_cori_term + metric_term + diffusion_info -> friction_acc[i];
        		}
            	if (h_index < NO_OF_VECTORS_V)
            	{
            		state_tendency -> velocity_gas[i] = forcings -> pot_vort_tend[i] - forcings -> e_kin_h_grad[i] + diffusion_info -> friction_acc[i];
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
            		state_tendency -> velocity_gas[i] = old_hor_grad_weight*interpolation -> pressure_gradient_acc_old[i] + new_hor_grad_weight*forcings -> pressure_gradient_acc[i] + forcings -> pot_vort_tend[i] - grid -> gravity_m[i] - forcings -> e_kin_h_grad[i] + hor_non_trad_cori_term + metric_term;
        		}
        		if (h_index < NO_OF_VECTORS_V)
        		{
            		state_tendency -> velocity_gas[i] = forcings -> pot_vort_tend[i] - forcings -> e_kin_h_grad[i];
        		}
        	}
        }
    }
    return 0;
}

int calc_partially_implicit_divvs(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
    scalar_times_vector(state_old -> density_dry, state_new -> velocity_gas, diagnostics -> mass_dry_flux_density, grid);
    divv_h(diagnostics -> mass_dry_flux_density, forcings -> mass_dry_flux_density_divv, grid);
    if (config_info -> scalar_diffusion_on == 1)
    {
        grad(state_old -> density_dry, diffusion_info -> mass_dry_diffusion_flux_density, grid);
        calc_mass_diffusion_coeffs(diagnostics -> temperature, state_old -> density_dry, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
        scalar_times_vector_scalar_h_v(diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v, diffusion_info -> mass_dry_diffusion_flux_density, diffusion_info -> mass_dry_diffusion_flux_density, grid);
        divv(diffusion_info -> mass_dry_diffusion_flux_density, diffusion_info -> mass_dry_diffusion_source_rate, grid, 0);
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    state_tendency -> density_dry[i] = -forcings -> mass_dry_flux_density_divv[i] + diffusion_info -> mass_dry_diffusion_source_rate[i];
	    }
    }
    else
    {
    	for (int i = 0; i < NO_OF_SCALARS; ++i)
        	state_tendency -> density_dry[i] = -forcings -> mass_dry_flux_density_divv[i];
    }
    scalar_times_vector(state_old -> t_tilde, state_new -> velocity_gas, diagnostics -> t_tilde_flux_density, grid);
    divv_h(diagnostics -> t_tilde_flux_density, forcings -> t_tilde_flux_density_divv, grid);
	divv_h(state_new -> velocity_gas, diagnostics -> wind_field_divv_h, grid);
    if (config_info -> rad_update == 1)
    {
        calc_rad_heating(radiation_tendency, NO_OF_SCALARS);
    }
    double rho_h, c_h_v;
    if (config_info -> scalar_diffusion_on == 1)
    {
        calc_temp_diffusion_coeffs(diagnostics -> temperature, state_old -> density_dry, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
        scalar_times_vector_scalar_h_v(diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v, diagnostics -> temp_gradient, diagnostics -> temperature_flux_density, grid);
        divv(diagnostics -> temperature_flux_density, diffusion_info -> temp_diffusion_heating, grid, 0);
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (config_info -> tracers_on == 1)
			{
				rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
				c_h_v = spec_heat_cap_diagnostics_v(state_old -> density_dry[i], state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
				diffusion_info -> temp_diffusion_heating[i] = rho_h*c_h_v*diffusion_info -> temp_diffusion_heating[i];
			}
			else
				diffusion_info -> temp_diffusion_heating[i] = state_old -> density_dry[i]*C_D_V*diffusion_info -> temp_diffusion_heating[i];
		}
    }

    double total_density;
    double c_v_cond;
    int h_index, layer_index;
    if (config_info -> tracers_on == 1)
    {
    	/*
    	phase transitions are on only at the third RK step
    	only then, they are also updated
    	*/
    	if (config_info -> phase_transitions_on == 1)
    	{
		    calc_h2otracers_source_rates(diffusion_info -> tracer_mass_source_rates, diffusion_info -> tracer_heat_source_rates, state_old-> tracer_densities, state_old -> tracer_density_temperatures, diagnostics -> temperature, NO_OF_TRACERS, NO_OF_SCALARS, delta_t);
		}
		for (int i = 0; i < NO_OF_TRACERS; ++i)
		{
		    for (int j = 0; j < NO_OF_SCALARS; ++j)
		        diffusion_info -> tracer_density[j] = state_old -> tracer_densities[i*NO_OF_SCALARS + j];
		    if (i < NO_OF_CONDENSATED_TRACERS)
		    {
		        for (int j = 0; j < NO_OF_VECTORS; ++j)
		        {
		            layer_index = j/NO_OF_VECTORS_PER_LAYER;
		            h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
		            if (h_index < NO_OF_VECTORS_V)
		                diffusion_info -> tracer_velocity[j] = state_new -> velocity_gas[j] - ret_sink_velocity(i, 0, 0.001);
		            else
		            	diffusion_info -> tracer_velocity[j] = state_new -> velocity_gas[j];
		        }
		        scalar_times_vector(diffusion_info -> tracer_density, diffusion_info -> tracer_velocity, diffusion_info -> tracer_flux_density, grid);
		        divv_h(diffusion_info -> tracer_flux_density, diffusion_info -> tracer_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
					diffusion_info -> tracer_density_temperature[j] = state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j];
		        scalar_times_vector(diffusion_info -> tracer_density_temperature, diffusion_info -> tracer_velocity, diffusion_info -> tracer_temperature_flux_density, grid);
		        divv_h(diffusion_info -> tracer_temperature_flux_density, diffusion_info -> tracer_temperature_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
				    c_v_cond = ret_c_v_cond(i, 0, state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j]/state_old -> tracer_densities[i*NO_OF_SCALARS + j]);
				    total_density = state_old -> density_dry[j];
				    for (int k = 0; k < NO_OF_TRACERS; ++k)
				        total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + j];
			        state_tendency -> tracer_density_temperatures[i*NO_OF_SCALARS + j] = -diffusion_info -> tracer_temperature_flux_density_divv[j] + state_old -> tracer_densities[i*NO_OF_SCALARS + j]/(c_v_cond*total_density)*(diffusion_info -> temp_diffusion_heating[j] + diffusion_info -> heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[i*NO_OF_SCALARS + j] + diffusion_info -> tracer_density_temperature[j]*config_info -> phase_transitions_on*(diffusion_info -> tracer_mass_source_rates[i*NO_OF_SCALARS + j]);
				}
		    }
		    else
		    {
		        scalar_times_vector_vector_h_v(diffusion_info -> tracer_density, state_new -> velocity_gas, state_new -> velocity_gas, diffusion_info -> tracer_flux_density, grid);
		        divv_h(diffusion_info -> tracer_flux_density, diffusion_info -> tracer_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
				    state_tendency -> tracer_densities[i*NO_OF_SCALARS + j] = -diffusion_info -> tracer_flux_density_divv[j] + config_info -> phase_transitions_on*diffusion_info -> tracer_mass_source_rates[i*NO_OF_SCALARS + j];
	            }
		    }
		}
    }
    double R_h, density_d_micro_value, density_v_micro_value, condensates_density_sum;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        if (config_info -> scalar_diffusion_on == 1)
        {
            total_density = state_old -> density_dry[i];
            for (int k = 0; k < NO_OF_TRACERS; ++k)
                total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + i];
            rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];		        
			condensates_density_sum = calc_condensates_density_sum(i, state_old -> tracer_densities);
			density_d_micro_value = calc_micro_density(state_old -> density_dry[i], condensates_density_sum);
			density_v_micro_value = calc_micro_density(state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i], condensates_density_sum);
		    c_h_v = spec_heat_cap_diagnostics_v(density_d_micro_value, density_v_micro_value);
		    R_h = gas_constant_diagnostics(density_d_micro_value, density_v_micro_value);
            state_tendency -> t_tilde[i] = -forcings -> t_tilde_flux_density_divv[i] - R_h/c_h_v*interpolation -> t_tilde[i]*diagnostics -> wind_field_divv_h[i] + 1/c_h_v*rho_h/total_density*(diffusion_info -> temp_diffusion_heating[i] + config_info -> momentum_diffusion_on*diffusion_info -> heating_diss[i] + radiation_tendency[i]) + config_info -> tracers_on*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
        }
        else
        {
            state_tendency -> t_tilde[i] = -forcings -> t_tilde_flux_density_divv[i] - R_D/C_D_V*interpolation -> t_tilde[i]*diagnostics -> wind_field_divv_h[i] + radiation_tendency[i]/C_D_V;
        }
    }
	scalar_times_vector(state_old -> entropy_gas, state_new -> velocity_gas, diagnostics -> entropy_gas_flux_density, grid);
	divv_h(diagnostics -> entropy_gas_flux_density, forcings -> entropy_gas_flux_density_divv, grid);
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
	    if (config_info -> scalar_diffusion_on == 1)
	    {
	        total_density = state_old -> density_dry[i];
	        for (int k = 0; k < NO_OF_TRACERS; ++k)
	            total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + i];
	        c_h_v = spec_heat_cap_diagnostics_v(state_old -> density_dry[i], state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
	        rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
	        state_tendency -> entropy_gas[i] = -forcings -> entropy_gas_flux_density_divv[i] + 1/diagnostics -> temperature[i]*(rho_h/total_density*(diffusion_info -> temp_diffusion_heating[i] + config_info -> momentum_diffusion_on*diffusion_info -> heating_diss[i] + radiation_tendency[i]) + config_info -> tracers_on*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
	    }
	    else
	    {
	        state_tendency -> entropy_gas[i] = -forcings -> entropy_gas_flux_density_divv[i] + 1/diagnostics -> temperature[i]*radiation_tendency[i];
	    }
	}
    return 0;
}


















