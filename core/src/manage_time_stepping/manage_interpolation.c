/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "../diagnostics/diagnostics.h"
#include <omp.h>
#include "../spatial_operators/spatial_operators.h"

int set_interpolated_temperature(State *state_old, State *state_new, Interpolate_info *interpolate, int totally_first_bool)
{
	double old_weight = R_D/C_D_P - 0.5;
	double new_weight = 1 - old_weight;
	if (totally_first_bool == 1)
	{
		old_weight = 1;
		new_weight = 0;
	}
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		interpolate -> temp_interpolate[i] = old_weight*state_old -> temp_gas[i] + new_weight*state_new -> temp_gas[i];
	}
	return 0;
}

int manage_pressure_gradient(State *current_state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	double old_hor_grad_weight = R_D/C_D_P - 0.5;
	double new_hor_grad_weight = 1 - old_hor_grad_weight;
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
	}
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
	// Before the calculation of the new pressure gradient, the old value needs to be stored for extrapolation.
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
	// The new pressure gradient os only calculated at the first RK step.
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		if (config_info -> tracers_on == 1)
			diagnostics -> c_h_p_field[i] = spec_heat_cap_diagnostics_p(current_state -> density_dry[i], current_state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
		else
			diagnostics -> c_h_p_field[i] = C_D_P;
	}
	scalar_times_grad(diagnostics -> c_h_p_field, current_state -> temp_gas, diagnostics -> pressure_gradient_0_m, grid);
	scalar_times_grad(interpolation -> temp_interpolate, diagnostics -> specific_entropy, interpolation -> pressure_gradient_1_interpolate, grid);
	// This is necessary for the vertical momentum equation.
	scalar_times_grad(current_state -> temp_gas, diagnostics -> specific_entropy, diagnostics -> pressure_gradient_1, grid);
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		forcings -> pressure_gradient_acc[i] = old_hor_grad_weight*(-interpolation -> pressure_gradient_0_old_m[i]) + new_hor_grad_weight*(-diagnostics -> pressure_gradient_0_m[i]) + interpolation -> pressure_gradient_1_interpolate[i];
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
	return 0;
}















	
