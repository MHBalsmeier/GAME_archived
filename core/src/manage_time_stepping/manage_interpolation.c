/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "../diagnostics/diagnostics.h"
#include <omp.h>
#include "../spatial_operators/spatial_operators.h"

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	
	// The weights for the horizontal pressure gradient acceleration extrapolation.
	double old_hor_grad_weight = -C_D_V/C_D_P;
	double new_hor_grad_weight = 1 - old_hor_grad_weight;
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
	}
	
	// Preparation of diagnostic quantities.
    double rho_h;
    #pragma omp parallel for private(rho_h)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	// Determining the speicific entropied of the dry air as well as of the water vapour.
		diagnostics -> specific_entropy_dry[i] = state -> entropy_density_dry[i]/(EPSILON_SECURITY + state -> density_dry[i]);
		diagnostics -> specific_entropy_vapour[i] = state -> tracer_entropy_densities[i]/(EPSILON_SECURITY + state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i]);
		// Determining the density of the gas phase.
		rho_h = state -> density_dry[i] + state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
		// The second pressure gradient terms for dry air as well as water vapour.
		diagnostics -> pressure_gradient_1_dry_prefactor[i] = state -> temperature_gas[i]*state -> density_dry[i]/rho_h;
		diagnostics -> pressure_gradient_1_vapour_prefactor[i] = state -> temperature_gas[i]*state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i]/rho_h;
	}
	
	// Before the calculation of the new pressure gradient, the old value needs to be stored for extrapolation.
	if (config_info -> totally_first_step_bool == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			interpolation -> pressure_gradient_0_old_m[i] = diagnostics -> pressure_gradient_0_m[i];
			interpolation -> pressure_gradient_1_old[i] = diagnostics -> pressure_gradient_1_dry[i] + diagnostics -> pressure_gradient_1_vapour[i];
		}
	}
	else
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			interpolation -> pressure_gradient_0_old_m[i] = 0;
			interpolation -> pressure_gradient_1_old[i] = 0;
		}
	}
	
	// Diagnozing c_h_p.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		if (config_info -> tracers_on == 1)
		{
			diagnostics -> c_h_p_field[i] = spec_heat_cap_diagnostics_p(state -> density_dry[i], state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i]);
		}
		else
		{
			diagnostics -> c_h_p_field[i] = C_D_P;
		}
	}
	
	// Calculating the three constituents of the pressure gradient acceleration.
	scalar_times_grad(diagnostics -> c_h_p_field, state -> temperature_gas, diagnostics -> pressure_gradient_0_m, grid);
	scalar_times_grad(diagnostics -> pressure_gradient_1_dry_prefactor, diagnostics -> specific_entropy_dry, diagnostics -> pressure_gradient_1_dry, grid);
	scalar_times_grad(diagnostics -> pressure_gradient_1_vapour_prefactor, diagnostics -> specific_entropy_vapour, diagnostics -> pressure_gradient_1_vapour, grid);
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		forcings -> pressure_gradient_acc[i] = 
		old_hor_grad_weight*(-interpolation -> pressure_gradient_0_old_m[i] + interpolation -> pressure_gradient_1_old[i])
		+ new_hor_grad_weight*(-diagnostics -> pressure_gradient_0_m[i] + diagnostics -> pressure_gradient_1_dry[i] + diagnostics -> pressure_gradient_1_vapour[i]);
	}
	
	// The pressure gradient has to get a deceleration factor in presence of condensates.
	double total_density;
	if (config_info -> tracers_on == 1)
	{
		#pragma omp parallel for private(rho_h, total_density)
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			rho_h = state -> density_dry[i] + state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
		    total_density = state -> density_dry[i];
		    for (int k = 0; k < NO_OF_TRACERS; ++k)
		    {
		        total_density += state -> tracer_densities[k*NO_OF_SCALARS + i];
	        }
			diffusion_info -> pressure_gradient_decel_factor[i] = rho_h/total_density;
		}
		scalar_times_vector(diffusion_info -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc, forcings -> pressure_gradient_acc, grid, 0);
	}
	
	return 0;
}















	
