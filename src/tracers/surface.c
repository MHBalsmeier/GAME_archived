/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating everything related to the surface.
*/

#include <stdlib.h>
#include <math.h>
#include "tracers.h"

double psi(double, double);

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




