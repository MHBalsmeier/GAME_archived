/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../../enum_and_typedefs.h"
#include "../../settings.h"
#include <stdlib.h>
#include <stdio.h>
#include "../../diagnostics/diagnostics.h"
#include <omp.h>
#include "../../spatial_operators/spatial_operators.h"

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Diffusion_info *diffusion_info, Config_info *config_info, int no_step_rk)
{
	
	// The weights for the horizontal pressure gradient acceleration extrapolation.
	double old_hor_grad_weight = -spec_heat_capacities_v_gas(0)/spec_heat_capacities_p_gas(0);
	double new_hor_grad_weight = 1 - old_hor_grad_weight;
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
	}
	
	// Preparation of diagnostic quantities.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		diagnostics -> pressure_gradient_1[i] = 0;
	}
    for (int j = 0; j < NO_OF_GASEOUS_CONSTITUENTS; ++j)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			// Determining the speicific entropied of the dry air as well as of the water vapour.
			diagnostics -> specific_entropy[i] = state -> entropy_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i]/
			state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i];
			// The second pressure gradient term prefactors for dry air as well as water vapour.
			diagnostics -> pressure_gradient_1_prefactor[i] = state -> temperature_gas[i]
			*state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i]/density_total(state, i);
		}
		scalar_times_grad(diagnostics -> pressure_gradient_1_prefactor, diagnostics -> specific_entropy, diagnostics -> pressure_gradient_1_component, grid);
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> pressure_gradient_1[i] += diagnostics -> pressure_gradient_1_component[i];
		}
	}
	// Before the calculation of the new pressure gradient, the old value needs to be stored for extrapolation.
	if (config_info -> totally_first_step_bool == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			interpolation -> pressure_gradient_0_old_m[i] = diagnostics -> pressure_gradient_0_m[i];
			interpolation -> pressure_gradient_1_old[i] = diagnostics -> pressure_gradient_1[i];
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
	
	// Diagnozing c_g_p.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i);
	}
	
	// Calculating the three constituents of the pressure gradient acceleration.
	scalar_times_grad(diagnostics -> c_g_p_field, state -> temperature_gas, diagnostics -> pressure_gradient_0_m, grid);
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		forcings -> pressure_gradient_acc[i] = 
		old_hor_grad_weight*(-interpolation -> pressure_gradient_0_old_m[i] + interpolation -> pressure_gradient_1_old[i])
		+ new_hor_grad_weight*(-diagnostics -> pressure_gradient_0_m[i] + diagnostics -> pressure_gradient_1[i]);
	}
	
	// The pressure gradient has to get a deceleration factor due to condensates.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diffusion_info -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
	}
	scalar_times_vector(diffusion_info -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc, forcings -> pressure_gradient_acc, grid, 0);
	
	return 0;
}















	
