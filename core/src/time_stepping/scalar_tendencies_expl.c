/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This is the horizontal (explicit) part of the constituent integration.
*/

#include "../enum_and_typedefs.h"
#include "atmostracers.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include "stdio.h"
#include "stdlib.h"

int scalar_tendencies_expl(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, int no_rk_step, Vector_field wind_advect_tracer)
{

	// declaring needed variables
    int h_index, layer_index;
    double c_v_cond;
    
    // determining the weights for the RK stepping
    double old_weight, new_weight, tracer_heating, density_gas_weight, density_total_weight;
    new_weight = 1;
    if (no_rk_step == 1)
    {
    	new_weight = 0.5;
    }
	old_weight = 1 - new_weight;
    
	/*
	phase transitions are only updated at the first RK step
	*/
	if (NO_OF_CONSTITUENTS == 4)
	{
	    calc_h2otracers_source_rates(
	    irrev -> constituent_mass_source_rates,
	    irrev -> constituent_heat_source_rates,
	    state -> mass_densities,
	    state -> condensed_density_temperatures,
	    state -> temperature_gas,
	    NO_OF_SCALARS,
	    delta_t,
	    config_info -> assume_lte);
	}
	
	// loop over all constituents
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		// Separating the mass density of the constituent at hand.
		#pragma omp parallel for
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
		    diagnostics -> scalar_field_placeholder[j] = state -> mass_densities[i*NO_OF_SCALARS + j];
	    }
        
        // This is the mass advection, which needs to be carried out for all constituents.
        // -------------------------------------------------------------------------------
		// calling mass advection according to user input
		if (config_info -> mass_advection_order == 2)
		{
			scalar_times_vector(diagnostics -> scalar_field_placeholder, wind_advect_tracer, diagnostics -> flux_density, grid);
		}
		if (config_info -> mass_advection_order == 3)
		{
			advection_3rd_order(diagnostics -> scalar_field_placeholder, state -> velocity_gas, wind_advect_tracer,
			diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid, no_rk_step, delta_t);
		}
        divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
		// adding the tendencies in all grid boxes
		#pragma omp parallel for private(layer_index, h_index)
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			layer_index = j/NO_OF_SCALARS_H;
			h_index = j - layer_index*NO_OF_SCALARS_H;
			if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
			{
				state_tendency -> mass_densities[i*NO_OF_SCALARS + j]
				= old_weight*state_tendency -> mass_densities[i*NO_OF_SCALARS + j]
				+ new_weight*(
				// the advection
				-diagnostics -> flux_density_divv[j]
				// the phase transition rates
				+ irrev -> constituent_mass_source_rates[i*NO_OF_SCALARS + j]);
		    }
	    }
	    
		// explicit entropy integrations
		// -----------------------------
		if ((config_info -> assume_lte == 1 && i == NO_OF_CONDENSED_CONSTITUENTS)
		|| (config_info -> assume_lte == 0 && i >= NO_OF_CONDENSED_CONSTITUENTS))
		{
			// Determining the specific entropy of the constituent at hand.
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				if (state -> mass_densities[i*NO_OF_SCALARS + j] != 0)
				{
					diagnostics -> scalar_field_placeholder[j] =
					state -> entropy_densities[(i - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j]
					/state -> mass_densities[i*NO_OF_SCALARS + j];
				}
				else
				{
					diagnostics -> scalar_field_placeholder[j] = 0;
				}
			}
			// calling entropy advection according to user input
			if (config_info -> entropy_advection_order == 2)
			{
				scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
			}
			if (config_info -> mass_advection_order == 3)
			{
				advection_3rd_order(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, wind_advect_tracer,
				diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid, no_rk_step, delta_t);
			}
			divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for private(layer_index, h_index, tracer_heating, density_gas_weight, density_total_weight)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				layer_index = j/NO_OF_SCALARS_H;
				h_index = j - layer_index*NO_OF_SCALARS_H;
				if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
				{
					// determining the heating rate that comes from the tracers
					tracer_heating = 0;
					density_gas_weight = 0;
					density_total_weight = 0;
					if (config_info -> assume_lte == 0)
					{
						density_gas_weight = density_gas(state, j);
						density_total_weight = density_total(state, j);
						tracer_heating = irrev -> constituent_heat_source_rates[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j];
					}
					if (config_info -> assume_lte == 1)
					{
						density_gas_weight = state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j];
						density_total_weight = state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j];
						for (int k = 0; k < NO_OF_CONDENSED_CONSTITUENTS + 1; ++k)
						{
							tracer_heating += irrev -> constituent_heat_source_rates[k*NO_OF_SCALARS + j];
						}
					}
					state_tendency -> entropy_densities[(i - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j]
					= old_weight*state_tendency -> entropy_densities[(i - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j]
					+ new_weight*(
					// the advection (resolved transport)
					-diagnostics -> flux_density_divv[j]
					// the diabatic forcings
					// weighting factor
					+ state -> mass_densities[i*NO_OF_SCALARS + j]/density_total_weight*(
					// dissipation of molecular + turbulent momentum diffusion
					irrev -> heating_diss[j]
					// molecular + turbulent heat transport
					+ irrev -> temperature_diffusion_heating[j]
					// radiation
					+ radiation_tendency[j]
					// this has to be divided by the temperature (we ware in the entropy equation)
					)/state -> temperature_gas[j]
					// phase transitions
					+ tracer_heating*state -> mass_densities[i*NO_OF_SCALARS + j]/density_gas_weight
					/state -> temperature_gas[j]);
				 }
			}
		}
    
		// This is the integration of the "density x temperature" fields. It only needs to be done for condensed constituents.
		// -------------------------------------------------------------------------------------------------------------------
		if (i < NO_OF_CONDENSED_CONSTITUENTS && config_info -> assume_lte == 0)
		{
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> condensed_density_temperatures[i*NO_OF_SCALARS + j];
			}
			// The constituent velocity has already been calculated.
		    scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> velocity_gas, diagnostics -> flux_density, grid);
		    divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for private(layer_index, h_index, c_v_cond)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				layer_index = j/NO_OF_SCALARS_H;
				h_index = j - layer_index*NO_OF_SCALARS_H;
				if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
				{
					c_v_cond = ret_c_v_cond(i, 0, state -> condensed_density_temperatures[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + state -> mass_densities[i*NO_OF_SCALARS + j]));
					state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j]
					= old_weight*state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j]
					+ new_weight*(
					// the advection
					-diagnostics -> flux_density_divv[j]
					// the source terms
					+ state -> mass_densities[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + c_v_cond*density_total(state, j))
					*(irrev -> temperature_diffusion_heating[j] + irrev -> heating_diss[j] + radiation_tendency[j])
					+ 1/c_v_cond*irrev -> constituent_heat_source_rates[i*NO_OF_SCALARS + j]
					+ diagnostics -> scalar_field_placeholder[j]*(irrev -> constituent_mass_source_rates[i*NO_OF_SCALARS + j]));
				}
			}
		}
	} // constituent loop
	return 0;
}









