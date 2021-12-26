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

int soil_interaction(State *state, Soil *soil, Diagnostics *diagnostics, Forcings *forcings, Grid *grid, double delta_t)
{
	/*
	This function computes the interaction of the dynamical core with the soil.
	*/
	
	double delta_z_soil = grid -> z_t_const/NO_OF_SOIL_LAYERS;

	int soil_index;
	// loop over all horizontal cells over land
	#pragma omp parallel for private(soil_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		if (grid -> is_land[i] == 1)
		{
			double heat_flux_density[NO_OF_SOIL_LAYERS];
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS - 1; ++soil_layer_index)
			{
				heat_flux_density[soil_layer_index]
				= grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(soil -> temperature[i + (soil_layer_index + 1)*NO_OF_SCALARS_H]
				- soil -> temperature[i + soil_layer_index*NO_OF_SCALARS_H])/delta_z_soil;
			}
			heat_flux_density[NO_OF_SOIL_LAYERS - 1]
			= grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(grid -> t_const_soil
			- soil -> temperature[i + (NO_OF_SOIL_LAYERS - 1)*NO_OF_SCALARS_H])/delta_z_soil;
			
			// summing up the power densities and transforming into a temperature change
			soil -> temperature[i]
			// sensible heat flux
			+= (soil -> power_flux_density_sensible[i]
			// latent heat flux
			+ soil -> power_flux_density_latent[i]
			// shortwave inbound radiation
			+ forcings -> sfc_sw_in[i]
			// longwave outbound radiation
			- forcings -> sfc_lw_out[i]
			// heat conduction from below
			+ heat_flux_density[0])
			/(delta_z_soil*grid -> sfc_rho_c[i])*delta_t;
			
			// loop over all soil layers below the first layer
			for (int soil_layer_index = 1; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
			{
				// index of this soil grid point
				soil_index = i + soil_layer_index*NO_OF_SCALARS_H;
				
				soil -> temperature[soil_index]
				// heat conduction from above
				+= (-heat_flux_density[soil_layer_index - 1]				
				// heat conduction from below
				+ heat_flux_density[soil_layer_index])
				/(delta_z_soil*grid -> sfc_rho_c[i])*delta_t;
			}
		}
	}
	return 0;
}







