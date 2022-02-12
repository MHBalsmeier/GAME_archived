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

double roughness_length_from_u10(double);
double psi(double, double);

const double KARMAN = 0.4;

int tke_update(Irreversible_quantities *irrev, double delta_t, State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function updates the specific turbulent kinetic energy (TKE), unit: J/kg.
	*/
	
	double decay_constant;
	// loop over all scalar gridpoints
	#pragma omp parallel for private(decay_constant)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// decay constant, as derived from diffusion
		decay_constant = 8*pow(M_PI, 2)/grid -> mean_velocity_area*(irrev -> viscosity_div[i] + irrev -> viscosity_curl[i])/density_gas(state, i);
		
		// prognostic equation for TKE
		irrev -> tke[i] += delta_t*(
		// production through dissipation of resolved energy
		irrev -> heating_diss[i]/density_gas(state, i)
		// decay through molecular dissipation
		- decay_constant*irrev -> tke[i]);
		
		// clipping negative values
		if (irrev -> tke[i] < 0)
		{
			irrev -> tke[i] = 0;
		}
	}
	return 0;
}

int update_roughness_length(Grid *grid, Diagnostics *diagnostics)
{
	/*
	This function updated the roughness length over water.
	*/
	
	double u10, z_agl;
	#pragma omp parallel for private(u10, z_agl)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// only over the sea the roughness length is time-dependant (because of the waves)
		if (grid -> is_land[i] == 0)
		{
			// calculating the 10 m wind velocity from the logarithmic wind profile
			z_agl = grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i];
			u10 = pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + i], 0.5)
			*log(10/grid -> roughness_length[i])
			/log(z_agl/grid -> roughness_length[i]);
			
			// calculating the roughness length fom the wind velocity
			grid -> roughness_length[i] = roughness_length_from_u10(u10);
		}
	}
	
	return 0;
}

double vertical_viscosity(double tke)
{
	/*
	This function returns the vertical kinematic Eddy viscosity as a function of the specific TKE.
	*/
	
	double prop_constant = 0.4; // unit: m
	double result = prop_constant*pow(tke, 0.5);
	return result;
}

double roughness_length_from_u10(double u10)
{
	/*
	This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.
	*/
	
	// refer to Stensrud, Parameterization schemes (2007), p.130
	
	// empirically determined formula for the SWH
	double swh = 0.0248*pow(u10, 2);
	
	// empirically determined period of the waves
	double period = 0.729*u10;
	
	// deep-water gravity waves
	double wavelength = GRAVITY_MEAN_SFC_ABS*pow(period, 2)/(2*M_PI);
	
	// final result
	double roughness_length = 1200*swh*pow(swh/wavelength, 4.5);
	
	return roughness_length;
}

double scalar_flux_resistance(double wind_h_lowest_layer, double z_agl, double roughness_length)
{
	/*
	This function returns the surface flux resistance for scalar quantities.
	*/
	
	double result = 1.0/(KARMAN*(roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length) + EPSILON_SECURITY))*
	// neutral conditions
	(log(z_agl/roughness_length)
	// non-neutral conditions
	- psi(z_agl, 100)
	// interfacial sublayer
	+ log(7));
	
	// limitting the result for security
	if (result < 50.0)
	{
		result = 50.0;
	}
	
	return result;
}

double momentum_flux_resistance(double wind_h_lowest_layer, double z_agl, double roughness_length)
{
	/*
	This function returns the surface flux resistance for momentum.
	*/
	
	double result = 1.0/(KARMAN*(roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length) + EPSILON_SECURITY))*
	// neutral conditions
	(log(z_agl/roughness_length)
	// non-neutral conditions
	- psi(z_agl, 100));
	
	// limitting the result for security
	if (result < 50.0)
	{
		result = 50.0;
	}
	
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








