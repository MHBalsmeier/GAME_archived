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
double roughness_length_from_u10(double);
double psi(double, double);

const double KARMAN = 0.4;

int tke_update(Irreversible_quantities *irrev, double delta_t, State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function updates the specific turbulent kinetic energy (TKE), unit: J/kg.
	*/
	
	// some constants
	double boundary_layer_height = 2500.0; // height of the boundary layer
	double mean_roughness_length = 0.6; // approximate global mean of the roughness length
	double roughness_length_exp = 1.0/2; // exponent of the roughness length
	double turb_prefactor = 2; // coefficient modulating the strength of the turbulence in the boundary layer
	// the e-folding time of TKE approximation in the boundary layer
	double tke_approx_time = 10800*pow(4, 5 - RES_ID);
	
	// think carefully before you change anything below this point
	
	// global integrals over the TKE and KE above the boundary layer
	double tke_glob_int_free = 0;
	double ke_glob_int_free = 0;
	double m_glob_int_free = 0;
	int i;
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			if (grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS + h_index] > boundary_layer_height)
			{
				tke_glob_int_free += irrev -> tke[i]*state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*grid -> volume[i];
				ke_glob_int_free += 0.5*diagnostics -> v_squared[i]*state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*grid -> volume[i];
				m_glob_int_free += state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*grid -> volume[i];
			}
		}
	}
	tke_glob_int_free = tke_glob_int_free/m_glob_int_free;
	ke_glob_int_free = ke_glob_int_free/m_glob_int_free;
	
	// the ratio of global unresolved to resolved kinetic energy above the boundary layer
	double tke_ke_ratio = tke_glob_int_free/(ke_glob_int_free + EPSILON_SECURITY);
	tke_ke_ratio = tke_ke_ratio;
	
	// computing the advection
	grad(irrev -> tke, diagnostics -> vector_field_placeholder, grid);
	inner_product(diagnostics -> vector_field_placeholder, state -> wind, diagnostics -> scalar_field_placeholder, grid);
	
	// loop over all scalar gridpoints
	double decay_constant, production_rate, ke, tke_expect, tke_expect_prefactor, u10, z_agl;
	#pragma omp parallel for private(i, decay_constant, production_rate, ke, tke_expect, tke_expect_prefactor, u10, z_agl)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		// updating the roughness length over water
		if (grid -> is_land[h_index] == 0)
		{
			// calculating the 10 m wind velocity from the logarithmic wind profile
			z_agl = grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + h_index] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index];
			u10 = pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + h_index], 0.5)
			*log(10/grid -> roughness_length[h_index])
			/log(z_agl/grid -> roughness_length[h_index]);
			
			// calculating the roughness length fom the wind velocity
			grid -> roughness_length[h_index] = roughness_length_from_u10(u10);
		}
		
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			
			// decay constant, as derived from diffusion
			decay_constant = 8*pow(M_PI, 2)/grid -> mean_velocity_area*(irrev -> viscosity_div[i] + irrev -> viscosity_curl[i])/density_gas(state, i);
			
			// generation of TKE in the boundary layer
			production_rate = 0;
			z_agl = grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index];
			if (z_agl <= boundary_layer_height)
			{
				// kinetic energy in this gridbox
				ke = 0.5*diagnostics -> v_squared[i];
				
				// expected value for the TKE from the energy spectrum assuming no boundary layer effects
				tke_expect = tke_ke_ratio*ke;
				
				tke_expect_prefactor = 1
				// factor taking into account the roughness of the surface
				+ turb_prefactor*pow(grid -> roughness_length[h_index]/mean_roughness_length, roughness_length_exp)
				// height-dependent factor
				*(exp(-z_agl/boundary_layer_height) - exp(-1.0))/(1 - exp(-1.0));
				
				// the amount of TKE that can be expected in this gridbox
				tke_expect = tke_expect_prefactor*tke_expect;
				
				// finally calculating the production rate
				production_rate = (tke_expect - irrev -> tke[i])/tke_approx_time;
				
				// restricting the production rate to positive values
				production_rate = fmax(0, production_rate);
			}
			
			// prognostic equation for TKE
			irrev -> tke[i] += delta_t*(
			// advection
			-diagnostics -> scalar_field_placeholder[i]
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








