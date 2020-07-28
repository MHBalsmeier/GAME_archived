/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this source file, the calculation of the tendencies, i.e. the explicit part of the dycore is managed.
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "atmostracers.h"
#include <stdlib.h>
#include <stdio.h>
#include "../../diagnostics/diagnostics.h"

int forward_tendencies(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	double old_hor_grad_weight = R_D/C_D_P - 0.5;
	double new_hor_grad_weight = 1 - old_hor_grad_weight;
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
	}
    if (config_info -> momentum_diffusion_on == 1 && no_step_rk == 2)
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
    // Now, the generalized Coriolis term is evaluated.
    coriolis_gen(diagnostics -> mass_dry_flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid);
    // Horizontal kinetic energy is prepared for the gradient term of the Lamb transformation.
    kinetic_energy(current_state -> velocity_gas, diagnostics -> e_kin_h, grid, 0);
    grad(diagnostics -> e_kin_h, forcings -> e_kin_h_grad, grid);
    // Now the explicit forces are added up.
    double metric_term, vertical_velocity, hor_non_trad_cori_term;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (i < NO_OF_SCALARS_H || i >= NO_OF_VECTORS - NO_OF_SCALARS_H)
            state_tendency -> velocity_gas[i] = 0;
        else
        {
        	if (config_info -> momentum_diffusion_on == 1 && no_step_rk == 2)
        	{
        		if (h_index >= NO_OF_SCALARS_H)
        		{
        			recov_hor_ver_pri(current_state -> velocity_gas, layer_index, h_index - NO_OF_SCALARS_H, &vertical_velocity, grid);
        			metric_term = -vertical_velocity*current_state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
        			hor_non_trad_cori_term = -vertical_velocity*dualgrid -> f_vec[2*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
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
        			hor_non_trad_cori_term = -vertical_velocity*dualgrid -> f_vec[2*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
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
















