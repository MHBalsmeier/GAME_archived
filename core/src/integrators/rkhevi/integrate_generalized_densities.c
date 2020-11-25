/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

/*
This is the horizontal (explicit) part of the constituent integration.
*/

#include "../../enum_and_typedefs.h"
#include "atmostracers.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../../diagnostics/diagnostics.h"
#include "stdio.h"
#include "stdlib.h"

int integrate_generalized_densities(State *state, Interpolation_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_rk_step)
{

    int h_index, layer_index;
    double c_v_cond, density_gas_value;
    
	/*
	phase transitions are on only at the first RK step
	only then, they are also updated
	*/
	if (no_rk_step == 0 && NO_OF_CONSTITUENTS == 4)
	{
	    calc_h2otracers_source_rates(irreversible_quantities -> constituent_mass_source_rates, irreversible_quantities -> constituent_heat_source_rates, state -> mass_densities, state -> condensed_density_temperatures, state -> temperature_gas, NO_OF_SCALARS, delta_t);
	}
	
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		// Separating the density of the constituent at hand.
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
					
					// The solid case.
					if (i < NO_OF_SOLID_CONSTITUENTS)
					{
						diagnostics -> velocity_gen[j] -= ret_sink_velocity(0, 0.001, density_gas_value);
					}
					// The liquid case.
					else
					{
						diagnostics -> velocity_gen[j] -= ret_sink_velocity(1, 0.001, density_gas_value);
					}
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
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
		    state_tendency -> mass_densities[i*NO_OF_SCALARS + j] =
			// the advection
		    -diagnostics -> flux_density_divv[j]
		    // the phase transition rates
		    + irreversible_quantities -> constituent_mass_source_rates[i*NO_OF_SCALARS + j];
	    }
		// Explicit entropy integrations
		// -----------------------------
		// Determining the specific entropy of the constituent at hand.
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			if (state -> mass_densities[i*NO_OF_SCALARS + j] != 0)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> entropy_densities[i*NO_OF_SCALARS + j]/state -> mass_densities[i*NO_OF_SCALARS + j];
			}
			else
			{
				diagnostics -> scalar_field_placeholder[j] = 0;
			}
		}
	    scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
	    divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			state_tendency -> entropy_densities[i*NO_OF_SCALARS + j] = 
			// the advection
			-diagnostics -> flux_density_divv[j]
			// the heating rates
			 + state -> mass_densities[i*NO_OF_SCALARS + j]/density_total(state, j)*radiation_tendency[j]/state -> temperature_gas[j];
	    }
    
		// This is the integration of the "density x temperature" fields. It only needs to be done for condensed constituents.
		if (i < NO_OF_CONDENSED_CONSTITUENTS)
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> condensed_density_temperatures[i*NO_OF_SCALARS + j];
			}
			// The constituent velocity has already been calculated.
		    scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> velocity_gen, diagnostics -> flux_density, grid);
		    divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				c_v_cond = ret_c_v_cond(i, 0, state -> condensed_density_temperatures[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + state -> mass_densities[i*NO_OF_SCALARS + j]));
			    state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j] =
			    // the advection
			    -diagnostics -> flux_density_divv[j]
			    // the source terms
			    + state -> mass_densities[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + c_v_cond*density_total(state, j))*(irreversible_quantities -> temperature_diffusion_heating[j] + irreversible_quantities -> heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*irreversible_quantities -> constituent_heat_source_rates[i*NO_OF_SCALARS + j] + diagnostics -> scalar_field_placeholder[j]*(irreversible_quantities -> constituent_mass_source_rates[i*NO_OF_SCALARS + j]);
			}
		}
	}
	return 0;
}









