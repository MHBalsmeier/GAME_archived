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

int scalar_tendencies_expl(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, int no_rk_step)
{

	// declaring needed variables
    int h_index, layer_index, k;
    double c_v_cond, density_gas_value, latent_heating;
    
    // determining the weights for the RK stepping
    double old_weight, new_weight;
    new_weight = 1;
    if (no_rk_step == 1)
    {
    	new_weight = 0.5;
    }
	old_weight = 1 - new_weight;
    
	/*
	phase transitions are only updated at the first RK step
	*/
	if (no_rk_step == 0 && NO_OF_CONSTITUENTS == 4)
	{
	    calc_h2otracers_source_rates(
	    irrev -> constituent_mass_source_rates,
	    irrev -> constituent_heat_source_rates,
	    state -> mass_densities,
	    state -> condensed_density_temperatures,
	    state -> temperature_gas,
	    NO_OF_SCALARS,
	    delta_t);
	}
	
	// loop over all constituents
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		// Separating the density of the constituent at hand.
		#pragma omp parallel for
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
		    diagnostics -> scalar_field_placeholder[j] = state -> mass_densities[i*NO_OF_SCALARS + j];
	    }
        
        // This is the mass advection, which needs to be carried out for all constituents.
        // -------------------------------------------------------------------------------
        // For condensed constituents, a sink velocity must be added.
        if (i < NO_OF_CONDENSED_CONSTITUENTS)
        {
	    	// Adding a sink velocity.
			#pragma omp parallel for private(layer_index, h_index, density_gas_value)
	        for (int j = 0; j < NO_OF_VECTORS; ++j)
	        {
	            layer_index = j/NO_OF_VECTORS_PER_LAYER;
	            h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
	            
	            // If the component is vertical, the sink velocity of this constituent must substracted.
	            if (h_index < NO_OF_SCALARS_H)
	            {
					if (layer_index == 0)
					{
						density_gas_value = density_gas(state, h_index);
					}
					else if (layer_index == NO_OF_LAYERS)
					{
						density_gas_value = density_gas(state, (layer_index - 1)*NO_OF_SCALARS_H + h_index);
					}
					else
					{
						density_gas_value = 0.5*(density_gas(state, (layer_index - 1)*NO_OF_SCALARS_H + h_index) + density_gas(state, layer_index*NO_OF_SCALARS_H + h_index));
					}
					diagnostics -> velocity_gen[j] -= ret_sink_velocity(0, 0.001, density_gas_value);
				}
                // The horizontal movement is not affected by the sink velocity.
	            else
	            {
	            	diagnostics -> velocity_gen[j] = state -> velocity_gas[j];
            	}
	        }
        	scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> velocity_gen, diagnostics -> flux_density, grid);
    	}
    	// This is not the case for gaseous constituents.
    	else
    	{
        	scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> velocity_gas, diagnostics -> flux_density, grid);
    	}
        divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
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
	    
		// Explicit entropy integrations
		if ((config_info -> assume_lte == 1 && i == NO_OF_CONDENSED_CONSTITUENTS)
		|| (config_info -> assume_lte == 0 && i >= NO_OF_CONDENSED_CONSTITUENTS))
		{
			// -----------------------------
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
			scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
			divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			#pragma omp parallel for private(layer_index, h_index, latent_heating, k)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				layer_index = j/NO_OF_SCALARS_H;
				h_index = j - layer_index*NO_OF_SCALARS_H;
				if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
				{
					// determing the latent heating depending on the configuration
					latent_heating = 0;
					if (config_info -> assume_lte == 0)
					{
						latent_heating = irrev -> constituent_heat_source_rates[i*NO_OF_SCALARS + j];
					}
					if (config_info -> assume_lte == 1)
					{
						// in this case, all the latent heating rates are assumed to act onto the gas phase
						for (k = 0; k < NO_OF_CONDENSED_CONSTITUENTS + 1; ++k)
						{
							latent_heating += irrev -> constituent_heat_source_rates[k*NO_OF_SCALARS + j];
						}
					}
					state_tendency -> entropy_densities[(i - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j]
					= old_weight*state_tendency -> entropy_densities[(i - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j]
					+ new_weight*(
					// the advection (resolved transport)
					-diagnostics -> flux_density_divv[j]
					// the diabatic forcings
					// weighting factor
					+ state -> mass_densities[i*NO_OF_SCALARS + j]/density_total(state, j)
					*(
					// dissipation of molecular + turbulent momentum diffusion
					irrev -> heating_diss[j]
					// molecular + turbulent heat transport
					+ irrev -> temperature_diffusion_heating[j]
					// radiation
					+ radiation_tendency[j]
					// phase transitions
					+ latent_heating
					// this has to be divided by the temperature (we ware in the entropy equation)
					)/state -> temperature_gas[j]);
				 }
			}
		}
    
		// This is the integration of the "density x temperature" fields. It only needs to be done for condensed constituents.
		if (i < NO_OF_CONDENSED_CONSTITUENTS && config_info -> assume_lte == 0)
		{
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> condensed_density_temperatures[i*NO_OF_SCALARS + j];
			}
			// The constituent velocity has already been calculated.
		    scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> velocity_gen, diagnostics -> flux_density, grid);
		    divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
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
	}
	return 0;
}









