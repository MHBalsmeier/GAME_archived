/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

/*
This is the horizontal (explicit) part of the tracer integration.
*/

#include "../../enum_and_typedefs.h"
#include "atmostracers.h"
#include "../../spatial_operators/spatial_operators.h"
#include "stdio.h"

int integrate_tracers(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
    double total_density;
    double c_v_cond;
    int h_index, layer_index;
	/*
	phase transitions are on only at the third RK step
	only then, they are also updated
	*/
	if (config_info -> phase_transitions_on == 1)
	{
	    calc_h2otracers_source_rates(diffusion_info -> tracer_mass_source_rates, diffusion_info -> tracer_heat_source_rates, state_old-> tracer_densities, state_old -> tracer_density_temperatures, state_old -> temperature_gas, NO_OF_TRACERS, NO_OF_SCALARS, delta_t);
	}
	
	// Determining the mass density of the tracer at hand.
	for (int i = 0; i < NO_OF_TRACERS; ++i)
	{
		if (no_rk_step == 0)
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
			    diffusion_info -> tracer_density[j] = state_old -> tracer_densities[i*NO_OF_SCALARS + j];
		    }
        }
		else
		{
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
			    diffusion_info -> tracer_density[j] = state_new -> tracer_densities[i*NO_OF_SCALARS + j];
		    }
        }
        
        // This is the mass advection, which needs to be carried out for all tracers.
        // For condensed tracers, a sink velocity must be added.
        if (i < NO_OF_CONDENSED_TRACERS)
        {
	    	// Adding a sink velocity.
	        for (int j = 0; j < NO_OF_VECTORS; ++j)
	        {
	            layer_index = j/NO_OF_VECTORS_PER_LAYER;
	            h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
	            // If the component is vertical, the sink velocity of this tracer must substracted.
	            if (h_index < NO_OF_SCALARS_H)
	            {
	                diffusion_info -> tracer_velocity[j] = state_new -> velocity_gas[j] - ret_sink_velocity(i, 0, 0.001);
                }
                // The horizontal movement is not affected by the sink velocity.
	            else
	            {
	            	diffusion_info -> tracer_velocity[j] = state_new -> velocity_gas[j];
            	}
	        }
        	scalar_times_vector_vector_h_v(diffusion_info -> tracer_density, diffusion_info -> tracer_velocity, diffusion_info -> tracer_velocity, diffusion_info -> tracer_flux_density, grid);
    	}
    	// This is not the case for gaseous tracers.
    	else
    	{
        	scalar_times_vector_vector_h_v(diffusion_info -> tracer_density, state_new -> velocity_gas, state_new -> velocity_gas, diffusion_info -> tracer_flux_density, grid);
    	}
        divv_h(diffusion_info -> tracer_flux_density, diffusion_info -> tracer_flux_density_divv, grid);
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
		    state_tendency -> tracer_densities[i*NO_OF_SCALARS + j] = 
		    // the advection
		    -diffusion_info -> tracer_flux_density_divv[j]
		    // the phase transition rates
		    + config_info -> phase_transitions_on*diffusion_info -> tracer_mass_source_rates[i*NO_OF_SCALARS + j];
		    // limiter
		    if (state_old -> tracer_densities[i*NO_OF_SCALARS + j] + delta_t*state_tendency -> tracer_densities[i*NO_OF_SCALARS + j] < 0)
		    {
		    	state_tendency -> tracer_densities[i*NO_OF_SCALARS + j] = -state_old -> tracer_densities[i*NO_OF_SCALARS + j]/delta_t;
		    }
        }
        
        // Explicit entropy integration (only needs to be done for the gaseous tracers)
        // Determining the entropy density of the tracer at hand.
        if (i >= NO_OF_CONDENSED_TRACERS)
        {
			if (no_rk_step == 0)
			{
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
					diffusion_info -> tracer_entropy_density[j] = state_old -> tracer_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j];
				}
			}
			else
			{
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
					diffusion_info -> tracer_entropy_density[j] = state_new -> tracer_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j];
				}
			}
		    scalar_times_vector_vector_h_v(diffusion_info -> tracer_entropy_density, state_new -> velocity_gas, state_new -> velocity_gas, diffusion_info -> tracer_flux_density, grid);
		    divv_h(diffusion_info -> tracer_flux_density, diffusion_info -> tracer_flux_density_divv, grid);
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				state_tendency -> tracer_entropy_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j] = 
				// the advection
				-diffusion_info -> tracer_flux_density_divv[j];
				// limiter
				if (state_old -> tracer_entropy_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j] + delta_t*state_tendency -> tracer_entropy_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j] < 0)
				{
					state_tendency -> tracer_entropy_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j] = -state_old -> tracer_entropy_densities[(i - NO_OF_CONDENSED_TRACERS)*NO_OF_SCALARS + j]/delta_t;
				}
		    }
        }
        
        // This is the integration of the "density x temperature" fields. It only needs to be done for condensed tracers.
	    if (i < NO_OF_CONDENSED_TRACERS)
	    {
	        
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diffusion_info -> tracer_density_temperature[j] = state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j];
			}
			// The tracer velocity has already been calculated.
	        scalar_times_vector(diffusion_info -> tracer_density_temperature, diffusion_info -> tracer_velocity, diffusion_info -> tracer_temperature_flux_density, grid, 0);
	        divv_h(diffusion_info -> tracer_temperature_flux_density, diffusion_info -> tracer_temperature_flux_density_divv, grid);
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
			    c_v_cond = ret_c_v_cond(i, 0, state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j]/state_old -> tracer_densities[i*NO_OF_SCALARS + j]);
			    total_density = state_old -> density_dry[j];
			    for (int k = 0; k < NO_OF_TRACERS; ++k)
			    {
			        total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + j];
		        }
		        state_tendency -> tracer_density_temperatures[i*NO_OF_SCALARS + j] =
		        // the advection
		         -diffusion_info -> tracer_temperature_flux_density_divv[j]
		        // the source terms
		         + state_old -> tracer_densities[i*NO_OF_SCALARS + j]/(c_v_cond*total_density)*(diffusion_info -> temperature_diffusion_heating[j] + diffusion_info -> heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[i*NO_OF_SCALARS + j] + diffusion_info -> tracer_density_temperature[j]*config_info -> phase_transitions_on*(diffusion_info -> tracer_mass_source_rates[i*NO_OF_SCALARS + j]);
		         
		      	// limiter
				if (state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j] + delta_t*state_tendency -> tracer_density_temperatures[i*NO_OF_SCALARS + j] < 0)
				{
					state_tendency -> tracer_density_temperatures[i*NO_OF_SCALARS + j] = -state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j]/delta_t;
				}
			}
	    }
	}
	
	return 0;
}









