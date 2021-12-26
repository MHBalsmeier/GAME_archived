/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains the soil component of GAME.
*/

#include <stdio.h>
#include <atmostracers.h>
#include "game_types.h"
#include "thermodynamics.h"
#include "time_stepping/time_stepping.h"

int soil_interaction(State *state, Soil *soil, Diagnostics *diagnostics, Forcings *forcings, Grid *grid, double delta_t)
{
	/*
	This function computes the interaction of the dynamical core with the soil.
	*/

	// The z-axis is positive upwards, negative downwards (as usual).	
	
	// loop over all horizontal cells over land
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// filtering for land points
		if (grid -> is_land[i] == 1)
		{
			// vectors needed for the Thomas algorithm
			double c_vector[NO_OF_SOIL_LAYERS - 1];
			double d_vector[NO_OF_SOIL_LAYERS];
			double e_vector[NO_OF_SOIL_LAYERS - 1];
			double r_vector[NO_OF_SOIL_LAYERS];
			double solution_vector[NO_OF_SOIL_LAYERS];
					
			// calculating the explicit part of the heat flux density
			double heat_flux_density_expl[NO_OF_SOIL_LAYERS];
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS - 1; ++soil_layer_index)
			{
				heat_flux_density_expl[soil_layer_index]
				= -grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(soil -> temperature[i + soil_layer_index*NO_OF_SCALARS_H]
				- soil -> temperature[i + (soil_layer_index + 1)*NO_OF_SCALARS_H])
				/(0.5*(grid -> z_soil_interface[soil_layer_index] - grid -> z_soil_interface[soil_layer_index + 2]));
			}
			heat_flux_density_expl[NO_OF_SOIL_LAYERS - 1]
			= -grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(soil -> temperature[i + (NO_OF_SOIL_LAYERS - 1)*NO_OF_SCALARS_H]
			- grid -> t_const_soil)
			/((grid -> z_soil_interface[NO_OF_SOIL_LAYERS - 1] + grid -> z_soil_interface[NO_OF_SOIL_LAYERS]) - 2*grid -> z_t_const);
			
			// calculating the explicit part of the temperature change
			r_vector[0]
			// old temperature
			= soil -> temperature[i]
			// sensible heat flux
			+ (soil -> power_flux_density_sensible[i]
			// latent heat flux
			+ soil -> power_flux_density_latent[i]
			// shortwave inbound radiation
			+ forcings -> sfc_sw_in[i]
			// longwave outbound radiation
			- forcings -> sfc_lw_out[i]
			// heat conduction from below
			+ heat_flux_density_expl[0])
			/((grid -> z_soil_interface[0] - grid -> z_soil_interface[1])*grid -> sfc_rho_c[i])*delta_t;
			
			// loop over all soil layers below the first layer
			for (int soil_layer_index = 1; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
			{
				
				r_vector[soil_layer_index]
				// old temperature
				= soil -> temperature[i + soil_layer_index*NO_OF_SCALARS_H]
				// heat conduction from above
				+ (-heat_flux_density_expl[soil_layer_index - 1]				
				// heat conduction from below
				+ heat_flux_density_expl[soil_layer_index])
				/((grid -> z_soil_interface[soil_layer_index] - grid -> z_soil_interface[soil_layer_index + 1])*grid -> sfc_rho_c[i])*delta_t;
			}
			
			// the diagonal component
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
			{
				d_vector[soil_layer_index] = 1.0;
			}
			// the off-diagonal components
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS - 1; ++soil_layer_index)
			{
				c_vector[soil_layer_index] = 0;
				e_vector[soil_layer_index] = 0;
			}
			
			// finally calling the Thomas algorithm to solve the system of linear equations
			thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_SOIL_LAYERS);
			
			// writing the result into the new state
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
			{
				soil -> temperature[i + soil_layer_index*NO_OF_SCALARS_H] = solution_vector[soil_layer_index];
			}
		
		}
	}
	return 0;
}







