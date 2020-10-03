/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "atmostracers.h"
#include "../../spatial_operators/spatial_operators.h"

int integrate_tracers(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
    double total_density;
    double c_v_cond;
    int h_index, layer_index;
    if (config_info -> tracers_on == 1)
    {
    	/*
    	phase transitions are on only at the third RK step
    	only then, they are also updated
    	*/
    	if (config_info -> phase_transitions_on == 1)
    	{
		    calc_h2otracers_source_rates(diffusion_info -> tracer_mass_source_rates, diffusion_info -> tracer_heat_source_rates, state_old-> tracer_densities, state_old -> tracer_density_temperatures, state_old -> temperature_gas, NO_OF_TRACERS, NO_OF_SCALARS, delta_t);
		}
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
		    if (i < NO_OF_CONDENSATED_TRACERS)
		    {
		        for (int j = 0; j < NO_OF_VECTORS; ++j)
		        {
		            layer_index = j/NO_OF_VECTORS_PER_LAYER;
		            h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
		            if (h_index < NO_OF_SCALARS_H)
		                diffusion_info -> tracer_velocity[j] = state_new -> velocity_gas[j] - ret_sink_velocity(i, 0, 0.001);
		            else
		            	diffusion_info -> tracer_velocity[j] = state_new -> velocity_gas[j];
		        }
		        scalar_times_vector(diffusion_info -> tracer_density, diffusion_info -> tracer_velocity, diffusion_info -> tracer_flux_density, grid);
		        divv_h(diffusion_info -> tracer_flux_density, diffusion_info -> tracer_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
					diffusion_info -> tracer_density_temperature[j] = state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j];
		        scalar_times_vector(diffusion_info -> tracer_density_temperature, diffusion_info -> tracer_velocity, diffusion_info -> tracer_temperature_flux_density, grid);
		        divv_h(diffusion_info -> tracer_temperature_flux_density, diffusion_info -> tracer_temperature_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
				    c_v_cond = ret_c_v_cond(i, 0, state_old -> tracer_density_temperatures[i*NO_OF_SCALARS + j]/state_old -> tracer_densities[i*NO_OF_SCALARS + j]);
				    total_density = state_old -> density_dry[j];
				    for (int k = 0; k < NO_OF_TRACERS; ++k)
				        total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + j];
			        state_tendency -> tracer_density_temperatures[i*NO_OF_SCALARS + j] = -diffusion_info -> tracer_temperature_flux_density_divv[j] + state_old -> tracer_densities[i*NO_OF_SCALARS + j]/(c_v_cond*total_density)*(diffusion_info -> temp_diffusion_heating[j] + diffusion_info -> heating_diss[j] + radiation_tendency[j]) + 1/c_v_cond*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[i*NO_OF_SCALARS + j] + diffusion_info -> tracer_density_temperature[j]*config_info -> phase_transitions_on*(diffusion_info -> tracer_mass_source_rates[i*NO_OF_SCALARS + j]);
				}
		    }
		    else
		    {
		        scalar_times_vector_vector_h_v(diffusion_info -> tracer_density, state_new -> velocity_gas, state_new -> velocity_gas, diffusion_info -> tracer_flux_density, grid);
		        divv_h(diffusion_info -> tracer_flux_density, diffusion_info -> tracer_flux_density_divv, grid);
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
				    state_tendency -> tracer_densities[i*NO_OF_SCALARS + j] = -diffusion_info -> tracer_flux_density_divv[j] + config_info -> phase_transitions_on*diffusion_info -> tracer_mass_source_rates[i*NO_OF_SCALARS + j];
	            }
		    }
		}
    }
	return 0;
}
