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

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Extrapolation_info *extrapolation, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_rk_step)
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
			extrapolation -> cpgradt_m_old[i] = 0;
		}
	}
	// This is the standard extrapolation.
	else
	{
		hor_pgrad_sound_extrapolation = get_impl_thermo_weight();
		old_hor_pgrad_sound_weight = -hor_pgrad_sound_extrapolation;
		new_hor_pgrad_sound_weight = 1 - old_hor_pgrad_sound_weight;
	}
	
	// The -c_p*grad(T) component is only updated at the zeroth RK step.
	if (no_rk_step == 0)
	{
		// Before calculating the new sound wave component of the pressure gradient acceleration, the old one must be saved for extrapolation.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			extrapolation -> cpgradt_m_old[i] =
			diagnostics -> cpgradt_m_cov_hor[i] + diagnostics -> cpgradt_cov_ver_corr_hor_m[i];
		}
		
		// 2.) The first pressure gradient term.
		// --------------------------------------------------------------------------------
		// Diagnozing c_g_p.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i);
		}
		// Multiplying c_g_p with the temperature gradient.
		grad_hor_cov(state -> temperature_gas, diagnostics -> cpgradt_m_cov_hor, grid);
		scalar_times_vector(diagnostics -> c_g_p_field, diagnostics -> cpgradt_m_cov_hor, diagnostics -> cpgradt_m_cov_hor, grid);
		
		// 3.) The first pressure gradient term. (vertical part and terrain correction)
		// --------------------------------------------------------------------------------
		// Diagnozing c_g_p.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i);
		}
		// Multplying c_g_p with the temperature gradient.
		grad_vert_cov(state -> temperature_gas, diagnostics -> cpgradt_cov_ver_corr_hor_m, grid);
		grad_oro_corr_no_add(diagnostics -> cpgradt_cov_ver_corr_hor_m, grid);
		scalar_times_vector(diagnostics -> c_g_p_field, diagnostics -> cpgradt_cov_ver_corr_hor_m, diagnostics -> cpgradt_cov_ver_corr_hor_m, grid);
	}
		
	// 4.) The second pressure gradient term. (updated at every step)
	// --------------------------------------------------------------------------------
	// Cleaning.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		diagnostics -> tgrads_cov_hor[i] = 0;
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
		grad_hor_cov(diagnostics -> scalar_field_placeholder, diagnostics -> pressure_gradient_1_component_cov, grid);
		scalar_times_vector(diagnostics -> pressure_gradient_1_prefactor, diagnostics -> pressure_gradient_1_component_cov, diagnostics -> pressure_gradient_1_component_cov, grid);
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> tgrads_cov_hor[i] += diagnostics -> pressure_gradient_1_component_cov[i];
		}
	}
		
	// 5.) The second pressure gradient term. (vertical part and terrain correction)
	// --------------------------------------------------------------------------------
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		diagnostics -> tgrads_cov_ver_corr_hor[i] = 0;
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
		grad_vert_cov(diagnostics -> scalar_field_placeholder, diagnostics -> pressure_gradient_1_component_corr, grid);
		grad_oro_corr_no_add(diagnostics -> pressure_gradient_1_component_corr, grid);
		scalar_times_vector(diagnostics -> pressure_gradient_1_prefactor, diagnostics -> pressure_gradient_1_component_corr, diagnostics -> pressure_gradient_1_component_corr, grid);
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> tgrads_cov_ver_corr_hor[i] += diagnostics -> pressure_gradient_1_component_corr[i];
		}
	}
	
    double expl_pgrad_sound_weight;
	expl_pgrad_sound_weight = 1 - get_impl_thermo_weight();
	// 6.) Here, the explicit part of the pressure gradient acceleration is added up.
	// --------------------------------------------------------------------------------
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
			old_hor_pgrad_sound_weight*(-extrapolation -> cpgradt_m_old[i])
			// the new time step horizontal -c_p*grad(T)
			+ new_hor_pgrad_sound_weight*(-diagnostics -> cpgradt_m_cov_hor[i] - diagnostics -> cpgradt_cov_ver_corr_hor_m[i])
			// covariant T*grad(s) component
			+ diagnostics -> tgrads_cov_hor[i]
			// the terrrain correction of the T*grad(s) component
			+ diagnostics -> tgrads_cov_ver_corr_hor[i];
		}
		// vertical case
		else
		{
			forcings -> pressure_gradient_acc_expl[i] =
			// weighted -c_p*grad(T) part
			expl_pgrad_sound_weight*(-diagnostics -> cpgradt_cov_ver_corr_hor_m[i])
			// fully explicit T*grad(s) part
			+ diagnostics -> tgrads_cov_ver_corr_hor[i];
		}
	}
	
	// 7.) The pressure gradient has to get a deceleration factor due to condensates.
	// --------------------------------------------------------------------------------
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irreversible_quantities -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
	}
	scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_expl, forcings -> pressure_gradient_acc_expl, grid);
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













	
