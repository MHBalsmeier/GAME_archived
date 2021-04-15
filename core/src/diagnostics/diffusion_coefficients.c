/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "diagnostics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int calc_mass_diffusion_coeffs(State *state, Config_info *config_info, Scalar_field mass_diffusion_coeff_numerical_h, Scalar_field mass_diffusion_coeff_numerical_v)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including Eddys).
	*/
	if (config_info -> mass_diff_h == 1 || config_info -> mass_diff_v == 1)
	{
		double mean_particle_mass = mean_particle_masses_gas(0);
		double eff_particle_radius = 130e-12;
		double mass_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v;
		#pragma omp parallel for private (mass_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v)
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &mass_diffusion_coeff);
		    if (config_info -> mass_diff_h == 1)
		    {
		    	upturn_for_scale_h = 1;
	    	}
		    else
		    {
		    	upturn_for_scale_h = 0;
	    	}
		    if (config_info -> mass_diff_v == 1)
		    {
		    	upturn_for_scale_v = 1;
	    	}
		    else
		    {
		    	upturn_for_scale_v = 0;
	    	}
		    mass_diffusion_coeff_numerical_h[i] = upturn_for_scale_h*mass_diffusion_coeff;
		    mass_diffusion_coeff_numerical_v[i] = upturn_for_scale_v*mass_diffusion_coeff;
		}
	}
	return 0;
}

int calc_temp_diffusion_coeffs(State *state, Config_info *config_info, Irreversible_quantities *irreversible_quantities, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including Eddys).
	*/
	// The Eddy viscosity coefficient only has to be calculated if it has not yet been done.
	if (config_info -> momentum_diff_h == 0)
	{
		hori_div_viscosity_eff(state, irreversible_quantities, grid, diagnostics, config_info, delta_t);
	}
	double c_g_v;
	#pragma omp parallel for private (c_g_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		c_g_v = spec_heat_cap_diagnostics_v(state, i, config_info);
		irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i] = c_g_v*irreversible_quantities -> viscosity_div_eff[i];
		// vertical Eddy viscosity is about two orders of magnitude smaller
		irreversible_quantities -> scalar_diffusion_coeff_numerical_v[i] = 0.01*irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i];
	}
	return 0;
}

int hori_div_viscosity_eff(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_cell*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	// the maximum diffusion coefficient (stability constraint)
	double max_diff_h_coeff_turb = 0.125*grid -> mean_area_cell/delta_t;
	double molecular_viscosity;
	
	#pragma omp parallel for private(molecular_viscosity)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// preliminary result
		irrev -> viscosity_div_eff[i] = 8*grid -> mean_area_cell*config_info -> diff_h_smag_fac
		*fabs(5.0/3*diagnostics -> velocity_gas_divv[i]);
		
		// calculating and adding the molecular viscosity
		calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i],
		eff_particle_radius, &molecular_viscosity);
		irrev -> viscosity_div_eff[i] += molecular_viscosity;
		
		// turbulent minimum
		if (irrev -> viscosity_div_eff[i] < min_diff_h_coeff_turb)
		{
			irrev -> viscosity_div_eff[i] = min_diff_h_coeff_turb;
		}
		
		// maximum (stability constraint)
		if (irrev -> viscosity_div_eff[i] > max_diff_h_coeff_turb)
		{
			irrev -> viscosity_div_eff[i] = max_diff_h_coeff_turb;
		}
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_div_eff[i] = density_gas(state, i)*irrev -> viscosity_div_eff[i];
	}
	return 0;
}

int hori_curl_viscosity_eff(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_cell*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	// the maximum diffusion coefficient (stability constraint)
	double max_diff_h_coeff_turb = 0.125*grid -> mean_area_cell/delta_t;
	double molecular_viscosity;
	
	int layer_index, h_index, scalar_index_from, scalar_index_to;
	#pragma omp parallel for private(molecular_viscosity, layer_index, h_index, scalar_index_from, scalar_index_to)
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		if (h_index >= NO_OF_SCALARS_H)
		{
			// preliminary result
			irrev -> viscosity_curl_eff[i] = 0.2*grid -> mean_area_cell*config_info -> diff_h_smag_fac
			*fabs(diagnostics -> rel_vort[NO_OF_VECTORS_H + 2*layer_index*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H]);
			
			// calculating and adding the molecular viscosity
			scalar_index_from = layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index - NO_OF_SCALARS_H];
			scalar_index_to = layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index - NO_OF_SCALARS_H];
			calc_diffusion_coeff(
			0.5*(state -> temperature_gas[scalar_index_from] + state -> temperature_gas[scalar_index_to]),
			mean_particle_mass,
			0.5*(state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_from] + state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_to]),
			eff_particle_radius, &molecular_viscosity);
			irrev -> viscosity_curl_eff[i] += molecular_viscosity;
			
			// turbulent minimum
			if (irrev -> viscosity_curl_eff[i] < min_diff_h_coeff_turb)
			{
				irrev -> viscosity_curl_eff[i] = min_diff_h_coeff_turb;
			}
			
			// maximum (stability constraint)
			if (irrev -> viscosity_curl_eff[i] > max_diff_h_coeff_turb)
			{
				irrev -> viscosity_curl_eff[i] = max_diff_h_coeff_turb;
			}
			
			// multiplying by the mass density of the gas phase
			irrev -> viscosity_curl_eff[i] = 0.5*(density_gas(state, scalar_index_from) + density_gas(state, scalar_index_to))*irrev -> viscosity_curl_eff[i];
		}
	}
	return 0;
}

int calc_diffusion_coeff(double temperature, double particle_mass, double density, double particle_radius, double *result)
{
    double thermal_velocity = sqrt(8*K_B*temperature/(M_PI*particle_mass));
    double particle_density = density/particle_mass;
    double cross_section = 4*M_PI*pow(particle_radius, 2);
    double mean_free_path = 1/(sqrt(2)*particle_density*cross_section);
    *result = 1.0/3.0*thermal_velocity*mean_free_path;
    return 0;
}







