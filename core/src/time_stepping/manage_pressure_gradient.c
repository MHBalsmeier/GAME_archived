/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include <stdlib.h>
#include <stdio.h>
#include "../diagnostics/diagnostics.h"
#include <omp.h>
#include "../spatial_operators/spatial_operators.h"

double pressure_gradient_1_damping_factor(double);

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Extrapolation_info *extrapolation, Irreversible_quantities *irreversible_quantities, Config_info *config_info)
{
	// 1.) The weights for the horizontal pressure gradient acceleration extrapolation.
	// --------------------------------------------------------------------------------
	double hor_pgrad_sound_extrapolation, old_hor_pgrad_sound_weight, new_hor_pgrad_sound_weight;
	// In the case of the first time step, no extrapolation is possible.
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_pgrad_sound_weight = 0;
		new_hor_pgrad_sound_weight = 1;
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			extrapolation -> pgrad_acc_old[i] = 0;
		}
	}
	// This is the standard extrapolation.
	else
	{
		hor_pgrad_sound_extrapolation = get_impl_thermo_weight();
		old_hor_pgrad_sound_weight = -hor_pgrad_sound_extrapolation;
		new_hor_pgrad_sound_weight = 1 - old_hor_pgrad_sound_weight;
	}
	
	// 2.) The first pressure gradient term (-c_p*grad(T)).
	// Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		extrapolation -> pgrad_acc_old[i] =	-diagnostics -> cpgradt[i] + diagnostics -> tgrads[i];
	}
	
	// --------------------------------------------------------------------------------
	// Diagnozing c_g_p.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i, config_info);
	}
	// Multiplying c_g_p with the temperature gradient.
	grad(state -> temperature_gas, diagnostics -> cpgradt, grid);
	scalar_times_vector(diagnostics -> c_g_p_field, diagnostics -> cpgradt, diagnostics -> cpgradt, grid);
		
	// 3.) The second pressure gradient term. (updated at every step)
	// --------------------------------------------------------------------------------
	// Cleaning.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		diagnostics -> tgrads[i] = 0;
	}
	// Each constitutent of the gas phase gets an individual term.
	int no_of_relevant_constituents;
	if (config_info -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
	}
	if (config_info -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
	}
	for (int j = 0; j < no_of_relevant_constituents; ++j)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			// Determining the speicific entropied of the dry air as well as of the water vapour.
			if (state -> mass_densities[(j + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i] != 0)
			{
				diagnostics -> scalar_field_placeholder[i] = state -> entropy_densities[j*NO_OF_SCALARS + i]/
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
		grad(diagnostics -> scalar_field_placeholder, diagnostics -> vector_field_placeholder, grid);
		scalar_times_vector(diagnostics -> pressure_gradient_1_prefactor, diagnostics -> vector_field_placeholder, diagnostics -> vector_field_placeholder, grid);
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> tgrads[i] += diagnostics -> vector_field_placeholder[i];
		}
	}
	
	// 4.) Here, the explicit part of the pressure gradient acceleration is added up.
	// --------------------------------------------------------------------------------
    double expl_pgrad_sound_weight;
	expl_pgrad_sound_weight = 1 - get_impl_thermo_weight();
	int layer_index, h_index;
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		// horizontal case
		if (h_index >= NO_OF_SCALARS_H)
		{
			forcings -> pressure_gradient_acc_expl[i] = 
			// the old time step horizontal -c_p*grad(T)
			old_hor_pgrad_sound_weight*extrapolation -> pgrad_acc_old[i]
			// the new time step horizontal -c_p*grad(T)
			+ new_hor_pgrad_sound_weight*(-diagnostics -> cpgradt[i] + diagnostics -> tgrads[i]);
		}
		// vertical case
		else
		{
			forcings -> pressure_gradient_acc_expl[i] = 
			// weighted -c_p*grad(T) part
			expl_pgrad_sound_weight*(-diagnostics -> cpgradt[i] + diagnostics -> tgrads[i]);
		}
	}
	
	// 5.) The pressure gradient has to get a deceleration factor due to condensates.
	// --------------------------------------------------------------------------------
	if (config_info -> assume_lte == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irreversible_quantities -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
		}
		scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_expl, forcings -> pressure_gradient_acc_expl, grid);
	}
	
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













	
