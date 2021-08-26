/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains the soil component of GAME.
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

// some for now constant parameters
const double thickness = 10;
// approximately the properties of water
const double density = 1000;
const double c_v = 4184;
const double heat_trans_coeff = 50;
const double t_min = 273.15 - 30;

int soil_interaction(Soil *soil, Diagnostics *diagnostics, Forcings *forcings, double delta_t)
{
	/*
	This function computes the interaction of the dynamical core with the soil.
	*/
	
	// loop over all horizontal cells
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// sensible heat flux density through the surface
		soil -> power_flux_density_sensible[i] = heat_trans_coeff*(diagnostics -> temperature_gas[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - soil -> temperature[i]);
		
		// summing up the power densities and transforming into a temperature change
		soil -> temperature[i]
		// sensible heat flux
		+= (soil -> power_flux_density_sensible[i]
		// shortwave inbound radiation
		+ forcings -> sfc_sw_in[i]
		// longwave outbound radiation
		- forcings -> sfc_lw_out[i])
		/(thickness*c_v*density)*delta_t;
		// clipping too low values
		if (soil -> temperature[i] < t_min)
		{
			soil -> temperature[i] = t_min;
		}
	}
		
	return 0;
}

int init_soil(Soil *soil, Diagnostics *diagnostics)
{
	/*
	This function initializes the soil state.
	*/
	
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// setting the soil temperature equal to the temperature in the lowest layer
		soil -> temperature[i] = diagnostics -> temperature_gas[NO_OF_SCALARS - NO_OF_SCALARS_H + i];
	}
	
	return 0;
}







