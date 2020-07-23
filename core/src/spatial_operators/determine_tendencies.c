/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this source file, the calculation of the tendencies, i.e. the explicit part of the dycore is managed.
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include "atmostracers.h"
#include "rte-rrtmgp-c.h"
#include <stdlib.h>
#include <stdio.h>

int integrate_continuity_dry(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_tracers(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_temp_gas(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_entropy_density_gas(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);

int explicit_momentum_tendencies(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	double old_hor_grad_weight = R_D/C_D_P - 0.5;
	double new_hor_grad_weight = 1 - old_hor_grad_weight;
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
	}
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
			diagnostics -> specific_entropy[i] = current_state -> entropy_density_gas[i]/rho_h;
    	}
    	else
    		diagnostics -> specific_entropy[i] = current_state -> entropy_density_gas[i]/current_state -> density_dry[i];
	}
    grad(diagnostics -> specific_entropy, diagnostics -> specific_entropy_gradient, grid);
	// Here, the pressure gradient extrapolation is managed.
	if (no_step_rk == 0)
	{
		if (config_info -> totally_first_step_bool == 0)
		{
			for (int i = 0; i < NO_OF_VECTORS; ++i)
			{
				interpolation -> pressure_gradient_0_old_m[i] = diagnostics -> pressure_gradient_0_m[i];
			}
		}
		else
		{
			for (int i = 0; i < NO_OF_VECTORS; ++i)
			{
				interpolation -> pressure_gradient_0_old_m[i] = 0;
			}
		}
	}
	scalar_times_vector(interpolation -> temp_interpolate, diagnostics -> specific_entropy_gradient, interpolation -> pressure_gradient_1_interpolate, grid);
	// This is necessary for the vertical momentum equation.
	scalar_times_vector(current_state -> temp_gas, diagnostics -> specific_entropy_gradient, diagnostics -> pressure_gradient_1, grid);
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	if (config_info -> tracers_on == 1)
    		diagnostics -> c_h_p_field[i] = spec_heat_cap_diagnostics_p(current_state -> density_dry[i], current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
    	else
    		diagnostics -> c_h_p_field[i] = C_D_P;
	}
    grad(current_state -> temp_gas, diagnostics -> temp_gradient, grid);
	scalar_times_vector(diagnostics -> c_h_p_field, diagnostics -> temp_gradient, diagnostics -> pressure_gradient_0_m, grid);
	if (no_step_rk == 0)
	{
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			forcings -> pressure_gradient_acc[i] = old_hor_grad_weight*(-interpolation -> pressure_gradient_0_old_m[i]) + new_hor_grad_weight*(-diagnostics -> pressure_gradient_0_m[i]) + interpolation -> pressure_gradient_1_interpolate[i];
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
        if (i < NO_OF_SCALARS_H || i >= NO_OF_VECTORS - NO_OF_SCALARS_H)
            state_tendency -> velocity_gas[i] = 0;
        else
        {
        	if (config_info -> momentum_diffusion_on == 1)
        	{
        		if (h_index >= NO_OF_SCALARS_H)
        		{
        			recov_hor_ver_pri(current_state -> velocity_gas, layer_index, h_index - NO_OF_SCALARS_H, &vertical_velocity, grid);
        			metric_term = -vertical_velocity*current_state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
        			hor_non_trad_cori_sign = -dualgrid -> h_curl_signs[4*(h_index - NO_OF_SCALARS_H) + 0];
        			hor_non_trad_cori_term = hor_non_trad_cori_sign*vertical_velocity*dualgrid -> f_vec[h_index - NO_OF_SCALARS_H];
            		state_tendency -> velocity_gas[i] = forcings -> pressure_gradient_acc[i] + forcings -> pot_vort_tend[i] - grid -> gravity_m[i] - forcings -> e_kin_h_grad[i] + hor_non_trad_cori_term + metric_term + diffusion_info -> friction_acc[i];
        		}
            	if (h_index < NO_OF_SCALARS_H)
            	{
            		state_tendency -> velocity_gas[i] = forcings -> pot_vort_tend[i] - forcings -> e_kin_h_grad[i] + diffusion_info -> friction_acc[i];
        		}
        	}
            else
            {
        		if (h_index >= NO_OF_SCALARS_H)
        		{
        			recov_hor_ver_pri(current_state -> velocity_gas, layer_index, h_index - NO_OF_SCALARS_H, &vertical_velocity, grid);
        			metric_term = -vertical_velocity*current_state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
        			hor_non_trad_cori_sign = -dualgrid -> h_curl_signs[4*(h_index - NO_OF_SCALARS_H) + 0];
        			hor_non_trad_cori_term = hor_non_trad_cori_sign*vertical_velocity*dualgrid -> f_vec[h_index - NO_OF_SCALARS_H];
            		state_tendency -> velocity_gas[i] = forcings -> pressure_gradient_acc[i] + forcings -> pot_vort_tend[i] - grid -> gravity_m[i] - forcings -> e_kin_h_grad[i] + hor_non_trad_cori_term + metric_term;
        		}
        		if (h_index < NO_OF_SCALARS_H)
        		{
            		state_tendency -> velocity_gas[i] = forcings -> pot_vort_tend[i] - forcings -> e_kin_h_grad[i];
        		}
        	}
        }
    }
    return 0;
}

int semi_implicit_scalar_tendencies(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	integrate_continuity_dry(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
    // Radiation is updated here.
    if (config_info -> rad_update == 1)
    {
        calc_rad_heating(radiation_tendency, NO_OF_SCALARS);
    }
    double rho_h, c_h_v;
    // Diffusion gets updated here.
    if (config_info -> scalar_diffusion_on == 1)
    {
        calc_temp_diffusion_coeffs(state_old -> temp_gas, state_old -> density_dry, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
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
	integrate_tracers(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
	integrate_temp_gas(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
	integrate_entropy_density_gas(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, diffusion_info, config_info, no_rk_step);
    return 0;
}

int integrate_continuity_dry(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
    scalar_times_vector(state_old -> density_dry, state_new -> velocity_gas, diagnostics -> mass_dry_flux_density, grid);
    divv_h(diagnostics -> mass_dry_flux_density, forcings -> mass_dry_flux_density_divv, grid);
    if (config_info -> scalar_diffusion_on == 1)
    {
        grad(state_old -> density_dry, diffusion_info -> mass_dry_diffusion_flux_density, grid);
        calc_mass_diffusion_coeffs(state_old -> temp_gas, state_old -> density_dry, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
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
	return 0;
}

int integrate_tracers(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
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
		    calc_h2otracers_source_rates(diffusion_info -> tracer_mass_source_rates, diffusion_info -> tracer_heat_source_rates, state_old-> tracer_densities, state_old -> tracer_density_temperatures, state_old -> temp_gas, NO_OF_TRACERS, NO_OF_SCALARS, delta_t);
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
		            if (h_index < NO_OF_SCALARS_H)
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
	return 0;
}

int integrate_temp_gas(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{   
    scalar_times_vector(state_old -> temp_gas, state_new -> velocity_gas, diagnostics -> temp_gas_flux, grid);
	divv_h(diagnostics -> temp_gas_flux, forcings -> temp_gas_flux_divv_h, grid);
	double total_density, rho_h;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        if (config_info -> scalar_diffusion_on == 1)
        {
            total_density = state_old -> density_dry[i];
            for (int k = 0; k < NO_OF_TRACERS; ++k)
                total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + i];
            rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
            state_tendency -> temp_gas[i] = -forcings -> temp_gas_flux_divv_h[i] + rho_h/total_density*(diffusion_info -> temp_diffusion_heating[i] + config_info -> momentum_diffusion_on*diffusion_info -> heating_diss[i] + radiation_tendency[i]) + config_info -> tracers_on*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
        }
        else
        {
            state_tendency -> temp_gas[i] = -forcings -> temp_gas_flux_divv_h[i] + radiation_tendency[i];
        }
    }
	return 0;
}

int integrate_entropy_density_gas(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	scalar_times_vector(state_old -> entropy_density_gas, state_new -> velocity_gas, diagnostics -> entropy_gas_flux_density, grid);
	divv_h(diagnostics -> entropy_gas_flux_density, forcings -> entropy_gas_flux_density_divv, grid);
	double rho_h, total_density;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
	    if (config_info -> scalar_diffusion_on == 1)
	    {
	        total_density = state_old -> density_dry[i];
	        for (int k = 0; k < NO_OF_TRACERS; ++k)
	            total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + i];
	        rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
	        state_tendency -> entropy_density_gas[i] = -forcings -> entropy_gas_flux_density_divv[i] + 1/state_old -> temp_gas[i]*(rho_h/total_density*(diffusion_info -> temp_diffusion_heating[i] + config_info -> momentum_diffusion_on*diffusion_info -> heating_diss[i] + radiation_tendency[i]) + config_info -> tracers_on*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
	    }
	    else
	    {
	        state_tendency -> entropy_density_gas[i] = -forcings -> entropy_gas_flux_density_divv[i] + 1/state_old -> temp_gas[i]*radiation_tendency[i];
	    }
	}
	return 0;
}
















