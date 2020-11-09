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

double pressure_gradient_1_damping_factor(double);

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolate_info *interpolation, Irreversible_quantities *irreversible_quantities, Config_info *config_info)
{
	// 1.) The weights for the horizontal pressure gradient acceleration extrapolation.
	// --------------------------------------------------------------------------------
	double old_hor_grad_weight, new_hor_grad_weight;
	// In the case of the first time step, no extrapolation is possible.
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_grad_weight = 0;
		new_hor_grad_weight = 1;
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			interpolation -> pressure_gradient_0_old_m[i] = 0;
			interpolation -> pressure_gradient_1_old[i] = 0;
		}
	}
	// This is the standard extrapolation.
	else
	{
		// Determining the weights.
		old_hor_grad_weight = -spec_heat_capacities_v_gas(0)/spec_heat_capacities_p_gas(0);
		new_hor_grad_weight = 1 - old_hor_grad_weight;
		// Before calculating the new pressure gradient acceleration, the old one must be saved for extrapolation.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			interpolation -> pressure_gradient_0_old_m[i] = diagnostics -> pressure_gradient_0_m[i];
			interpolation -> pressure_gradient_1_old[i] = diagnostics -> pressure_gradient_1[i];
		}
	}
	
	// 2.) The first pressure gradient term.
	// --------------------------------------------------------------------------------
	// Diagnozing c_g_p.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i);
	}
	// Multplying c_g_p with the temperature gradient.
	scalar_times_grad(diagnostics -> c_g_p_field, state -> temperature_gas, diagnostics -> pressure_gradient_0_m, grid);
	
	// 3.) The second pressure gradient term.
	// --------------------------------------------------------------------------------
	// Cleaning.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		diagnostics -> pressure_gradient_1[i] = 0;
	}
	// Each constitutent of the gas phase gets an individual term.
    for (int j = 0; j < NO_OF_GASEOUS_CONSTITUENTS; ++j)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			// Determining the speicific entropied of the dry air as well as of the water vapour.
			if (state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i] != 0)
			{
				diagnostics -> scalar_field_placeholder[i] = state -> entropy_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i]/
				state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i];
			}
			else
			{
				diagnostics -> scalar_field_placeholder[i] = 0;
			}
			// The second pressure gradient term prefactors for dry air as well as water vapour.
			diagnostics -> pressure_gradient_1_prefactor[i] =
			// damping term for small densities
			pressure_gradient_1_damping_factor(state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i])
			// the physical term
			*state -> temperature_gas[i]
			*state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i]/density_gas(state, i);
		}
		scalar_times_grad(diagnostics -> pressure_gradient_1_prefactor, diagnostics -> scalar_field_placeholder, diagnostics -> pressure_gradient_1_component, grid);
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> pressure_gradient_1[i] += diagnostics -> pressure_gradient_1_component[i];
		}
	}
	
	// 4.) Here, the pressure gradient acceleration is added up.
	// --------------------------------------------------------------------------------
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		forcings -> pressure_gradient_acc[i] = 
		old_hor_grad_weight*(-interpolation -> pressure_gradient_0_old_m[i] + interpolation -> pressure_gradient_1_old[i])
		+ new_hor_grad_weight*(-diagnostics -> pressure_gradient_0_m[i] + diagnostics -> pressure_gradient_1[i]);
	}
	
	// 5.) The pressure gradient has to get a deceleration factor due to condensates.
	// --------------------------------------------------------------------------------
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irreversible_quantities -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
	}
	scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc, forcings -> pressure_gradient_acc, grid, 0);
	
	return 0;
}

double pressure_gradient_1_damping_factor(double density_value)
{
	double safe_density = 1e-8;
	double result;
	result = density_value/safe_density;
	if (result > 1)
	{
		result = 1;
	}
	return result;
}













	
