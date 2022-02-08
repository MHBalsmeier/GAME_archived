/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, diffusion coefficients, including Eddy viscosities, are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"
#include "subgrid_scale.h"

double ver_hor_viscosity(double, double, double);
double swh_from_u10(double);
double roughness_length_from_swh(double);
double psi(double, double);

const double KARMAN = 0.4;

int tke_update(Irreversible_quantities *irrev, double delta_t, State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function updates the specific turbulent kinetic energy (TKE), unit: J/kg.
	*/
	// the ratio of global unresolved to resolved kinetic energy
	double tke_ke_ratio = 0.1*pow(4, 5 - RES_ID);
	// the e-folding time of TKE approximation
	double tke_approx_time = 10800*pow(4, 5 - RES_ID);
	// computing the advection
	grad(irrev -> tke, diagnostics -> vector_field_placeholder, grid);
	inner_product(diagnostics -> vector_field_placeholder, state -> wind, diagnostics -> scalar_field_placeholder, grid);
	double boundary_layer_height = 1000.0;
	double roughness_length_factor = 1.0/0.08;
	int i;
	double decay_constant, production_rate, ke, tke_expect, u10;
	#pragma omp parallel for private(i, decay_constant, production_rate, ke, tke_expect, u10)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		// updating the roughness length over water
		if (grid -> is_land[h_index] == 0)
		{
			u10 = pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + h_index], 0.5)
			*log(10/grid -> roughness_length[h_index])
			/log((grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + h_index] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index])/grid -> roughness_length[h_index]);
			grid -> roughness_length[h_index] = roughness_length_from_swh(swh_from_u10(u10));
		}
		
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			
			// decay constant, as derived from diffusion
			decay_constant = 8*pow(M_PI, 2)/grid -> mean_velocity_area*(irrev -> viscosity_div[i] + irrev -> viscosity_curl[i])/density_gas(state, i);
			
			// generation of TKE in the boundary layer
			production_rate = 0;
			if (grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index] <= boundary_layer_height)
			{
				// kinetic energy in this gridbox
				ke = 0.5*diagnostics -> v_squared[i];
				// expected value for the TKE from the energy spectrum
				tke_expect = tke_ke_ratio*ke;
				
				production_rate =
				// factor taking into account the roughness of the surface
				roughness_length_factor*grid -> roughness_length[h_index]
				// height-dependent factor
				*(boundary_layer_height - (grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index]))/boundary_layer_height
				*(tke_expect - irrev -> tke[i])/tke_approx_time;
				
				// restricting the production rate to positive values
				production_rate = fmax(0, production_rate);
			}
			
			// prognostic equation for TKE
			irrev -> tke[i] += delta_t*(
			// advection
			- diagnostics -> scalar_field_placeholder[i]
			// production through dissipation of resolved energy
			+ irrev -> heating_diss[i]/density_gas(state, i)
			// decay through molecular dissipation
			- decay_constant*irrev -> tke[i]
			// production through turbulence generation in the boundary layer
			+ production_rate);
			// clipping negative values which might occur through advection
			if (irrev -> tke[i] < 0)
			{
				irrev -> tke[i] = 0;
			}
		}
	}
	return 0;
}

double ver_hor_viscosity(double tke, double delta_z, double mixing_length)
{
	/*
	This function returns the vertical kinematic Eddy viscosity as a function of the specific TKE.
	*/
	
	double prop_constant = 0.004*fmin(delta_z, mixing_length); // unit: m
	double result = prop_constant*pow(tke, 0.5);
	return result;
}

double swh_from_u10(double u10)
{
	/*
	This function returns the significant wave height (SWH) as a function of the 10 m wind velocity.
	*/
	
	double swh_max = 12;
	
	// formula for the SWH
	double swh = u10/7;
	
	// restricting the maximum wave height
	if (swh > swh_max)
	{
		swh = swh_max;
	}
	
	return swh;
}


double roughness_length_from_swh(double swh)
{
	/*
	This function returns the roughness length as a function of the significant wave height (SWH).
	*/
	
	double roughness_length = swh/35;
	return roughness_length;
}

double sfc_flux_resistance(double wind_h_lowest_layer, double z_agl, double roughness_length)
{
	/*
	This function returns the surface flux resistance.
	*/
	
	double result = 1.0/(KARMAN*(roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length) + EPSILON_SECURITY))*
	// neutral conditions
	(log(z_agl/roughness_length)
	// non-neutral conditions
	- psi(z_agl, 100)
	// interfacial sublayer
	+ log(7));
	return result;
}

double roughness_velocity(double wind_speed, double z_agl, double roughness_length)
{
	/*
	This function returns the roughness velocity.
	*/
	
	double result = wind_speed*KARMAN/log(z_agl/roughness_length);
	return result;
}

double psi(double z_eff, double l)
{
	/*
	This is a helper function for the correction to the surface flux resistance for non-neutral conditions.
	*/
	
	// z_eff: effective height above the surface
	// l: Monin-Obukhov length
	
	// avoiding l == 0
	if (l == 0)
	{
		l = EPSILON_SECURITY;
	}
	
	// helper variable
	double x = pow(1 - 15*z_eff/l, 0.25);
	
	double result;
	// unstable conditions
	if (l < 0)
	{
		result = 2*log((1 + pow(x, 2))/2);
	}
	// neutral and stable conditions
	else
	{
		result = -4.7*z_eff/l;
	}
	return result;
}








