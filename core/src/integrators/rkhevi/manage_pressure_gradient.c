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

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Interpolation_info *interpolation, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_rk_step)
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
			interpolation -> pgrad_cpgradt_m_old[i] = 0;
			interpolation -> pgrad_tgrads_old[i] = 0;
		}
	}
	// This is the standard extrapolation.
	else
	{
		old_hor_grad_weight = -get_impl_thermo_weight();
		new_hor_grad_weight = 1 + get_impl_thermo_weight();
	}
	
	// The pressure gradient is only updated at the zeroth RK step.
	if (no_rk_step == 0)
	{
		// Before calculating the new pressure gradient acceleration, the old one must be saved for extrapolation.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			interpolation -> pgrad_cpgradt_m_old[i] =
			diagnostics -> pgrad_cpgradt_m_cov_hor[i] + diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m[i];
			interpolation -> pgrad_tgrads_old[i] =
			diagnostics -> pgrad_tgrads_cov_hor[i] + diagnostics -> pgrad_tgrads_cov_ver_corr_hor[i];
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
		grad_hor_cov(state -> temperature_gas, diagnostics -> pgrad_cpgradt_m_cov_hor, grid);
		scalar_times_vector(diagnostics -> c_g_p_field, diagnostics -> pgrad_cpgradt_m_cov_hor, diagnostics -> pgrad_cpgradt_m_cov_hor, grid);
		
		// 3.) The second pressure gradient term.
		// --------------------------------------------------------------------------------
		// Cleaning.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> pgrad_tgrads_cov_hor[i] = 0;
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
				diagnostics -> pgrad_tgrads_cov_hor[i] += diagnostics -> pressure_gradient_1_component_cov[i];
			}
		}
	
		// 4.) The first pressure gradient term. (vertical part and terrain correction)
		// --------------------------------------------------------------------------------
		// Diagnozing c_g_p.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			diagnostics -> c_g_p_field[i] = spec_heat_cap_diagnostics_p(state, i);
		}
		// Multplying c_g_p with the temperature gradient.
		grad_vert_cov(state -> temperature_gas, diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m, grid);
		grad_oro_corr_no_add(diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m, grid);
		scalar_times_vector(diagnostics -> c_g_p_field, diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m, diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m, grid);
			
		// 5.) The second pressure gradient term. (vertical part and terrain correction)
		// --------------------------------------------------------------------------------
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			diagnostics -> pgrad_tgrads_cov_ver_corr_hor[i] = 0;
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
				diagnostics -> pgrad_tgrads_cov_ver_corr_hor[i] += diagnostics -> pressure_gradient_1_component_corr[i];
			}
		}
		
		// 6.) Here, the pressure gradient acceleration is added up.
		// --------------------------------------------------------------------------------
		int layer_index, h_index;
		for (int i = 0; i < NO_OF_VECTORS; ++i)
		{
			layer_index = i/NO_OF_VECTORS_PER_LAYER;
			h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
			// horizontal case
			if (h_index >= NO_OF_SCALARS_H)
			{
				forcings -> pressure_gradient_acc[i] =
				// the old time step convariant gradient
				old_hor_grad_weight*(-interpolation -> pgrad_cpgradt_m_old[i] + interpolation -> pgrad_tgrads_old[i])
				// the new time step convariant gradient
				+ new_hor_grad_weight*(-diagnostics -> pgrad_cpgradt_m_cov_hor[i] + diagnostics -> pgrad_tgrads_cov_hor[i]
				// the terrrain correction
				- diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m[i] + diagnostics -> pgrad_tgrads_cov_ver_corr_hor[i]);
			}
			// vertical case
			else
			{
				forcings -> pressure_gradient_acc[i] =
				-diagnostics -> pgrad_cpgradt_cov_ver_corr_hor_m[i] + diagnostics -> pgrad_tgrads_cov_ver_corr_hor[i];
			}
		}
		
		// 7.) The pressure gradient has to get a deceleration factor due to condensates.
		// --------------------------------------------------------------------------------
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irreversible_quantities -> pressure_gradient_decel_factor[i] = density_gas(state, i)/density_total(state, i);
		}
		scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc, forcings -> pressure_gradient_acc, grid);
	}
	return 0;
}

double pressure_gradient_1_damping_factor(double density_value)
{
	double safe_density = 1e-6;
	double result;
	result = density_value/safe_density;
	if (result > 1)
	{
		result = 1;
	}
	return result;
}













	
