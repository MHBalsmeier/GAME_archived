/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
In this file, the explicit component of the pressure gradient acceleration is managed.
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include <stdlib.h>
#include <stdio.h>
#include "../thermodynamics/thermodynamics.h"
#include <omp.h>
#include "../spatial_operators/spatial_operators.h"

double pressure_gradient_1_damping_factor(double);

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Extrapolation_info *extrapolation, Irreversible_quantities *irrev, Config_info *config_info)
{
	/*
	This function computes the pressure gradient acceleration.
	*/
	
	// 1.) The weights for the horizontal pressure gradient acceleration extrapolation.
	// --------------------------------------------------------------------------------
	double hor_pgrad_sound_extrapolation, old_hor_pgrad_sound_weight, new_hor_pgrad_sound_weight;
	// In the case of the first time step, no extrapolation is possible.
	if (config_info -> totally_first_step_bool == 1)
	{
		old_hor_pgrad_sound_weight = 0;
		new_hor_pgrad_sound_weight = 1;
	}
	// This is the standard extrapolation.
	else
	{
		hor_pgrad_sound_extrapolation = get_impl_thermo_weight() - 0.5;
		old_hor_pgrad_sound_weight = -hor_pgrad_sound_extrapolation;
		new_hor_pgrad_sound_weight = 1 - old_hor_pgrad_sound_weight;
	}
	
	// 2.) the nonlinear pressure gradient term
	// Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		extrapolation -> pgrad_acc_old[i] =	forcings -> pressure_gradient_acc_nl_expl[i] + forcings -> pressure_gradient_acc_l_expl[i];
	}
	
	// diagnozing c_g_p and multiplying by the full potential tempertature
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i, config_info);
		diagnostics -> scalar_field_placeholder[i] = diagnostics -> c_g_p_field[i]*(grid -> theta_bg[i] + state -> theta_pert[i]);
	}
	// multiplying c_g_p by the temperature gradient
	grad(state -> exner_pert, forcings -> pressure_gradient_acc_nl_expl, grid);
	scalar_times_vector(diagnostics -> scalar_field_placeholder, forcings -> pressure_gradient_acc_nl_expl, forcings -> pressure_gradient_acc_nl_expl, grid);
		
	// 3.) the linear pressure gradient term
	// -------------------------------------
	// cleaning and diagnozing c_g_p*theta_pert
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		forcings -> pressure_gradient_acc_l_expl[i] = grid -> exner_bg_grad[i];
		diagnostics -> scalar_field_placeholder[i] = diagnostics -> c_g_p_field[i]*state -> theta_pert[i];
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, forcings -> pressure_gradient_acc_l_expl, forcings -> pressure_gradient_acc_l_expl, grid);
	
	
	// 4.) Here, the explicit part of the pressure gradient acceleration is added up.
	// --------------------------------------------------------------------------------
    double expl_pgrad_sound_weight;
	expl_pgrad_sound_weight = 1 - get_impl_thermo_weight();
	int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		// horizontal case
		if (h_index >= NO_OF_SCALARS_H)
		{
			forcings -> pressure_gradient_acc_expl[i] = old_hor_pgrad_sound_weight*extrapolation -> pgrad_acc_old[i]
			+ new_hor_pgrad_sound_weight*(-diagnostics -> cpgradt[i] + diagnostics -> tgrads[i]);
		}
		// vertical case
		else
		{
			forcings -> pressure_gradient_acc_expl[i] = expl_pgrad_sound_weight*(-diagnostics -> cpgradt[i] + diagnostics -> tgrads[i]);
		}
	}
	
	// 5.) The pressure gradient has to get a deceleration factor due to condensates.
	// --------------------------------------------------------------------------------
	if (config_info -> assume_lte == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irrev -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
		}
		scalar_times_vector(irrev -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_expl, forcings -> pressure_gradient_acc_expl, grid);
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













	
