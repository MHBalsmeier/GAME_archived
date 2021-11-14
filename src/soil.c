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

// one for now constant parameters
const double thickness = 1;

int soil_interaction(State *state, Soil *soil, Diagnostics *diagnostics, Forcings *forcings, Grid *grid, double delta_t)
{
	/*
	This function computes the interaction of the dynamical core with the soil.
	*/
	
	double flux_resistance;
	// loop over all horizontal cells
	#pragma omp parallel for private(flux_resistance)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// sensible heat flux density through the surface
		flux_resistance = sfc_flux_resistance(pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + i], 0.5), grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i]
		- grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i], grid -> roughness_length[i]);
		soil -> power_flux_density_sensible[i] = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + NO_OF_SCALARS - NO_OF_SCALARS_H + i]
		*spec_heat_capacities_v_gas(0)*(diagnostics -> temperature_gas[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - soil -> temperature[i])/flux_resistance;
		// summing up the power densities and transforming into a temperature change
		soil -> temperature[i]
		// sensible heat flux
		+= (soil -> power_flux_density_sensible[i]
		// shortwave inbound radiation
		+ forcings -> sfc_sw_in[i]
		// longwave outbound radiation
		- forcings -> sfc_lw_out[i])
		/(thickness*grid -> sfc_rho_c[i])*delta_t;
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







