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

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info)
{
	/*
	This function computes the pressure gradient acceleration.
	*/
	
	// 2.) the nonlinear pressure gradient term
	// Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
	if (config_info -> totally_first_step_bool == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			forcings -> pgrad_acc_old[i] = -forcings -> pressure_gradient_acc_neg_nl[i] - forcings -> pressure_gradient_acc_neg_l[i];
		}
	}
	
	// diagnozing c_g_p and multiplying by the full potential tempertature
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i, config_info);
		diagnostics -> scalar_field_placeholder[i] = diagnostics -> c_g_p_field[i]*(grid -> theta_bg[i] + state -> theta_pert[i]);
	}
	grad(state -> exner_pert, forcings -> pressure_gradient_acc_neg_nl, grid);
	scalar_times_vector(diagnostics -> scalar_field_placeholder, forcings -> pressure_gradient_acc_neg_nl, forcings -> pressure_gradient_acc_neg_nl, grid);
		
	// 3.) the linear pressure gradient term
	// -------------------------------------
	// cleaning and diagnozing c_g_p*theta_pert
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = diagnostics -> c_g_p_field[i]*state -> theta_pert[i];
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, grid -> exner_bg_grad, forcings -> pressure_gradient_acc_neg_l, grid);
	
	// 4.) The pressure gradient has to get a deceleration factor due to condensates.
	// --------------------------------------------------------------------------------
	if (config_info -> assume_lte == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irrev -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
		}
		scalar_times_vector(irrev -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_neg_nl, forcings -> pressure_gradient_acc_neg_nl, grid);
		scalar_times_vector(irrev -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_neg_l, forcings -> pressure_gradient_acc_neg_l, grid);
	}
	
	// at the very fist step, the old time step pressure gradient acceleration must be saved here
	if (config_info -> totally_first_step_bool == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			forcings -> pgrad_acc_old[i] = -forcings -> pressure_gradient_acc_neg_nl[i] - forcings -> pressure_gradient_acc_neg_l[i];
		}
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













	
