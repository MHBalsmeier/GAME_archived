/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, quantities referring to the planetary boundary layer are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "subgrid_scale.h"

double roughness_length_from_u10_sea(double);
double scalar_flux_resistance(double, double, double, double);
double psi_h(double, double);
double psi_m(double, double);

const double KARMAN = 0.4;

int update_sfc_turb_quantities(State *state, Grid *grid, Diagnostics *diagnostics, Config *config, double delta_t)
{
	/*
	This function updates surface-related turbulence quantities.
	*/
	
	double u_lowest_layer, u10, z_agl, theta_lowest_layer, theta_second_layer, dz, dtheta_dz, w_pert, theta_pert, w_pert_theta_pert_avg;
	// semi-empirical coefficient
	double prop_coeff = 0.2;
	#pragma omp parallel for private(u_lowest_layer, u10, z_agl, theta_lowest_layer, theta_second_layer, dz, dtheta_dz, w_pert, theta_pert, w_pert_theta_pert_avg)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		z_agl = grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i];
		
		// wind speed in the lowest layer
		u_lowest_layer = pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + i], 0.5);
			
		// calculating the 10 m wind velocity from the logarithmic wind profile
		u10 = u_lowest_layer*log(10.0/grid -> roughness_length[i])/log(z_agl/grid -> roughness_length[i]);
		
		// only over the sea the roughness length is time-dependant (because of the waves)
		if (grid -> is_land[i] == 0)
		{
			
			// calculating the roughness length fom the wind velocity
			grid -> roughness_length[i] = roughness_length_from_u10_sea(u10);
		}
		
		// updating the roughness velocity
		diagnostics -> roughness_velocity[i] = roughness_velocity(u_lowest_layer, z_agl, grid -> roughness_length[i]);
		
		// theta in the lowest layer
		theta_lowest_layer = grid  -> theta_bg[NO_OF_SCALARS - NO_OF_SCALARS_H + i] + state -> theta_pert[NO_OF_SCALARS - NO_OF_SCALARS_H + i];
		// theta in the second-lowest layer
		theta_second_layer = grid  -> theta_bg[NO_OF_SCALARS - 2*NO_OF_SCALARS_H + i] + state -> theta_pert[NO_OF_SCALARS - 2*NO_OF_SCALARS_H + i];
		
		// delta z
		dz = grid -> z_scalar[NO_OF_SCALARS - 2*NO_OF_SCALARS_H + i] - grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i];
		
		// vertical gradient of theta
		dtheta_dz = (theta_second_layer - theta_lowest_layer)/dz;
		
		// the perturbation of the vertical velocity is assumed to be proportional to the 10 m wind speed
		// times a stability-dependant factor
		w_pert = u10*fmax(0.001, 0.02*(1.0 - dtheta_dz/0.01));
		theta_pert = -0.2*delta_t*w_pert*dtheta_dz;
		w_pert_theta_pert_avg = prop_coeff*w_pert*theta_pert;
		
		// security
		if (fabs(w_pert_theta_pert_avg) < EPSILON_SECURITY)
		{
			w_pert_theta_pert_avg = EPSILON_SECURITY;
		}
		
		// the result
		diagnostics -> monin_obukhov_length[i] = -theta_lowest_layer*pow(diagnostics -> roughness_velocity[i], 3.0)/(KARMAN*G_MEAN_SFC_ABS*w_pert_theta_pert_avg);
	}
	
	// updating the surface flux resistance acting on scalar quantities (moisture and sensible heat)
	if (config -> soil_on == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			diagnostics -> scalar_flux_resistance[i] = scalar_flux_resistance(diagnostics -> roughness_velocity[i],
			grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i],
			grid -> roughness_length[i], diagnostics -> monin_obukhov_length[i]);
		}
	}
	
	return 0;
}

double roughness_length_from_u10_sea(double u10)
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
	double wavelength = G_MEAN_SFC_ABS*pow(period, 2)/(2*M_PI);
	
	// final result
	double roughness_length = 1200*swh*pow(swh/fmax(wavelength, EPSILON_SECURITY), 4.5);
	
	// avoid too small values for stability
	return fmax(0.001, roughness_length);
}

double scalar_flux_resistance(double roughness_velocity_value, double z_agl, double roughness_length_value, double monin_obukhov_length_value)
{
	/*
	This function returns the surface flux resistance for scalar quantities.
	*/
	
	double result = 1.0/(KARMAN*roughness_velocity_value)*
	// neutral conditions
	(log(z_agl/roughness_length_value)
	// non-neutral conditions
	- psi_h(z_agl, monin_obukhov_length_value)
	// interfacial sublayer
	+ log(7.0));
	
	// limitting the result for security
	if (result < 50.0)
	{
		result = 50.0;
	}
	
	return result;
}

double momentum_flux_resistance(double wind_h_lowest_layer, double z_agl, double roughness_length_value, double monin_obukhov_length_value)
{
	/*
	This function returns the surface flux resistance for momentum.
	*/
	
	double result = 1.0/(KARMAN*roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length_value))*
	// neutral conditions
	(log(z_agl/roughness_length_value)
	// non-neutral conditions
	- psi_m(z_agl, monin_obukhov_length_value));
	
	// limitting the result for security
	if (result < 50.0)
	{
		result = 50.0;
	}
	
	return result;
}

double roughness_velocity(double wind_speed, double z_agl, double roughness_length_value)
{
	/*
	This function returns the roughness velocity.
	*/
	
	double denominator = log(z_agl/roughness_length_value);
	
	// security
	if (fabs(denominator) < EPSILON_SECURITY)
	{
		denominator = EPSILON_SECURITY;
	}
	
	double result = wind_speed*KARMAN/denominator;
	
	return fmax(EPSILON_SECURITY, result);
}

double psi_h(double z_eff, double l)
{
	/*
	This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.
	*/
	
	// z_eff: effective height above the surface
	// l: Monin-Obukhov length
	
	// avoiding l == 0
	if (fabs(l) < EPSILON_SECURITY)
	{
		l = EPSILON_SECURITY;
	}
	
	double result;
	// unstable conditions
	if (l < 0.0)
	{
		// helper variable
		double x = pow(1.0 - 15.0*z_eff/l, 0.25);
		
		result = 2.0*log((1.0 + pow(x, 2.0))/2.0);
	}
	// neutral and stable conditions
	else
	{
		result = -4.7*z_eff/l;
	}
	
	return result;
}

double psi_m(double z_eff, double l)
{
	/*
	This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.
	*/
	
	// z_eff: effective height above the surface
	// l: Monin-Obukhov length
	
	// avoiding l == 0
	if (fabs(l) < EPSILON_SECURITY)
	{
		l = EPSILON_SECURITY;
	}
	
	double result;
	// unstable conditions
	if (l < 0.0)
	{
		// helper variable
		double x = pow(1.0 - 15.0*z_eff/l, 0.25);
		
		result = 2.0*log((1 + x)/2.0) + log((1.0 + pow(x, 2.0))/2.0) - 2.0*atan(x) + M_PI/2.0;
	}
	// neutral and stable conditions
	else
	{
		result = -4.7*z_eff/l;
	}
	
	return result;
}








