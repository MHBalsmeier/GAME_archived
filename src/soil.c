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

// one for now constant parameter
const double thickness = 1;

int soil_interaction(State *state, Soil *soil, Diagnostics *diagnostics, Forcings *forcings, Grid *grid, double delta_t)
{
	/*
	This function computes the interaction of the dynamical core with the soil.
	*/
	
	// loop over all horizontal cells over land
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		if (grid -> is_land[i] == 1)
		{
			// summing up the power densities and transforming into a temperature change
			soil -> temperature[i]
			// sensible heat flux
			+= (soil -> power_flux_density_sensible[i]
			// latent heat flux
			+ soil -> power_flux_density_latent[i]
			// shortwave inbound radiation
			+ forcings -> sfc_sw_in[i]
			// longwave outbound radiation
			- forcings -> sfc_lw_out[i])
			/(thickness*grid -> sfc_rho_c[i])*delta_t;
		}
	}
	return 0;
}







