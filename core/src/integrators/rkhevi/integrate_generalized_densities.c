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

int integrate_generalized_densities(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{

    int h_index, layer_index;
    double c_v_cond, density_gas_value;
    
	/*
	phase transitions are on only at the third RK step
	only then, they are also updated
	*/
	if (config_info -> phase_transitions_on == 1)
	{
	    calc_h2otracers_source_rates(diffusion_info -> constituent_mass_source_rates, diffusion_info -> constituent_heat_source_rates, state_old-> mass_densities, state_old -> condensed_density_temperatures, state_old -> temperature_gas, NO_OF_CONSTITUENTS, NO_OF_SCALARS, delta_t);
	}
	
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		// Separating the density of the constituent at hand.
		if (no_rk_step == 0)
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
			    diagnostics -> density_gen[j] = state_old -> mass_densities[i*NO_OF_SCALARS + j];
		    }
        }
		else
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
			    diagnostics -> density_gen[j] = state_new -> mass_densities[i*NO_OF_SCALARS + j];
		    }
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
						density_gas_value = density_gas(state_old, h_index);
					}
					else if (layer_index == NO_OF_LAYERS)
					{
						density_gas_value = density_gas(state_old, (layer_index - 1)*NO_OF_SCALARS + h_index);
					}
					else
					{
						density_gas_value = 0.5*(density_gas(state_old, (layer_index - 1)*NO_OF_SCALARS + h_index) + density_gas(state_old, layer_index*NO_OF_SCALARS + h_index));
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
	            	diagnostics -> velocity_gen[j] = state_new -> velocity_gas[j];
            	}
	        }
        	scalar_times_vector_vector_h_v(diagnostics -> density_gen, diagnostics -> velocity_gen, diagnostics -> velocity_gen, diagnostics -> flux_density, grid);
    	}
    	// This is not the case for gaseous constituents.
    	else
    	{
        	scalar_times_vector_vector_h_v(diagnostics -> density_gen, state_new -> velocity_gas, state_new -> velocity_gas, diagnostics -> flux_density, grid);
    	}
        divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
		    state_tendency -> mass_densities[i*NO_OF_SCALARS + j] = 
		    // the advection
		    -diagnostics -> flux_density_divv[j]
		    // the phase transition rates
		    + config_info -> phase_transitions_on*diffusion_info -> constituent_mass_source_rates[i*NO_OF_SCALARS + j];
		    // limiter
		    if (state_old -> mass_densities[i*NO_OF_SCALARS + j] + delta_t*state_tendency -> mass_densities[i*NO_OF_SCALARS + j] < 0)
		    {
		    	state_tendency -> mass_densities[i*NO_OF_SCALARS + j] = -state_old -> mass_densities[i*NO_OF_SCALARS + j]/delta_t;
		    }
        }
        
		// Explicit entropy integrations
		// -----------------------------
		// Determining the entropy density of the constituent at hand.
		if (no_rk_step == 0)
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> density_gen[j] = state_old -> entropy_densities[i*NO_OF_SCALARS + j];
			}
		}
		else
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> density_gen[j] = state_new -> entropy_densities[i*NO_OF_SCALARS + j];
			}
		}
	    scalar_times_vector_vector_h_v(diagnostics -> density_gen, state_new -> velocity_gas, state_new -> velocity_gas, diagnostics -> flux_density, grid);
	    divv_h(diagnostics -> flux_density, diagnostics-> flux_density_divv, grid);
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			state_tendency -> entropy_densities[i*NO_OF_SCALARS + j] = 
			// the advection
			-diagnostics -> flux_density_divv[j];
			// limiter
			if (state_old -> entropy_densities[i*NO_OF_SCALARS + j] + delta_t*state_tendency -> entropy_densities[i*NO_OF_SCALARS + j] < 0)
			{
				state_tendency -> entropy_densities[i*NO_OF_SCALARS + j] = -state_old -> entropy_densities[i*NO_OF_SCALARS + j]/delta_t;
			}
	    }
    
		// This is the integration of the "density x temperature" fields. It only needs to be done for condensed constituents.
		if (i < NO_OF_CONDENSED_CONSTITUENTS)
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> density_gen[j] = state_old -> condensed_density_temperatures[i*NO_OF_SCALARS + j];
			}
			// The constituent velocity has already been calculated.
		    scalar_times_vector(diagnostics -> density_gen, diagnostics -> velocity_gen, diagnostics -> flux_density, grid, 0);
		    divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				c_v_cond = ret_c_v_cond(i, 0, state_old -> condensed_density_temperatures[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + state_old -> mass_densities[i*NO_OF_SCALARS + j]));
			    state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j] =
			    // the advection
			     -diagnostics -> flux_density_divv[j]
			    // the source terms
			    + state_old -> mass_densities[i*NO_OF_SCALARS + j]/(c_v_cond*density_total(state_old, j))*(diffusion_info -> temperature_diffusion_heating[j] + diffusion_info -> heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*config_info -> phase_transitions_on*diffusion_info -> constituent_heat_source_rates[i*NO_OF_SCALARS + j] + diagnostics -> density_gen[j]*config_info -> phase_transitions_on*(diffusion_info -> constituent_mass_source_rates[i*NO_OF_SCALARS + j]);
			     
			  	// limiter
				if (state_old -> condensed_density_temperatures[i*NO_OF_SCALARS + j] + delta_t*state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j] < 0)
				{
					state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j] = -state_old -> condensed_density_temperatures[i*NO_OF_SCALARS + j]/delta_t;
				}
			}
		}
	}
	
	return 0;
}









