/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "diagnostics.h"
#include <stdlib.h>
#include <stdio.h>


int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int calc_mass_diffusion_coeffs(State *state, Config_info *config_info, Scalar_field mass_diffusion_coeff_numerical_h, Scalar_field mass_diffusion_coeff_numerical_v)
{
	if (config_info -> mass_diff_h == 1 || config_info -> mass_diff_v == 1)
	{
		double mean_particle_mass = mean_particle_masses_gas(0);
		double eff_particle_radius = 130e-12;
		double mass_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v;
		#pragma omp parallel for private (mass_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v)
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &mass_diffusion_coeff);
			// homogeneous for now
			calc_diffusion_coeff(273.15, mean_particle_mass, 1, eff_particle_radius, &mass_diffusion_coeff);
		    if (config_info -> mass_diff_h == 1)
		    {
		    	upturn_for_scale_h = 1.5*pow(10, 4);
	    	}
		    else
		    {
		    	upturn_for_scale_h = 0;
	    	}
		    if (config_info -> mass_diff_v == 1)
		    {
		    	upturn_for_scale_v = 1.5*pow(10, 3);
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

int calc_temp_diffusion_coeffs(State *state, Config_info *config_info, Scalar_field temp_diffusion_coeff_numerical_h, Scalar_field temp_diffusion_coeff_numerical_v)
{
	double mean_particle_mass = mean_particle_masses_gas(0);
	double eff_particle_radius = 130e-12;
	double temp_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v, rho_g, c_g_v;
	#pragma omp parallel for private (temp_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v, rho_g, c_g_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
	    calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &temp_diffusion_coeff);
		// homogeneous for now
		calc_diffusion_coeff(273.15, mean_particle_mass, 1, eff_particle_radius, &temp_diffusion_coeff);
	    if (config_info -> temperature_diff_h == 1)
	    {
			upturn_for_scale_h = 1.5*pow(10, 4);
    	}
		else
	    {
			upturn_for_scale_h = 0;
    	}
	    if (config_info -> temperature_diff_v == 1)
	    {
			upturn_for_scale_v = 1.5*pow(10, 3);
    	}
		else
	    {
			upturn_for_scale_v = 0;
    	}
		rho_g = density_gas(state, i);
		c_g_v = spec_heat_cap_diagnostics_v(state, i, config_info);
	    temp_diffusion_coeff_numerical_h[i] = upturn_for_scale_h*rho_g*c_g_v*temp_diffusion_coeff;
	    temp_diffusion_coeff_numerical_v[i] = upturn_for_scale_v*rho_g*c_g_v*temp_diffusion_coeff;
	}
	return 0;
}

int calc_divv_term_viscosity_eff(State *state, Scalar_field divv_term_viscosity_eff, Grid *grid, double delta_t)
{
	// these things can be modified
	double diff_h_smag_fac = 0.18;
	double shear_bg = 5*1.5e-5;
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	double molecular_viscosity;
	// the maximum diffusion coefficient (based on stability, including a 50 % safety margin)
	double max_diff_h_coeff_turb = 0.5*0.25*grid -> mean_area_edge/delta_t;
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_edge*diff_h_smag_fac*shear_bg;
	double viscosity_value;
	#pragma omp parallel for private(molecular_viscosity)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// calculating the molecular viscosity
		// calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &molecular_viscosity);
		// homogeneous for now
		viscosity_value = min_diff_h_coeff_turb;
		calc_diffusion_coeff(273.15, mean_particle_mass, 1, eff_particle_radius, &molecular_viscosity);
		// the molecular viscosity is the absolute minimum
		if (molecular_viscosity > viscosity_value)
		{
			viscosity_value = molecular_viscosity;
		}
		// neglecting the volume viscosity for now
		divv_term_viscosity_eff[i] = 4.0/3*viscosity_value;
	}
	return 0;
}

int calc_curl_term_viscosity_eff(State *state, Scalar_field curl_term_viscosity_eff, Grid *grid, double delta_t)
{
	// these things can be modified
	double diff_h_smag_fac = 0.18;
	double shear_bg = 5*1.5e-5;
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	double molecular_viscosity;
	// the maximum diffusion coefficient (based on stability, including a 50 % safety margin)
	double max_diff_h_coeff_turb = 0.5*0.25*grid -> mean_area_edge/delta_t;
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_edge*diff_h_smag_fac*shear_bg;
	double viscosity_value;
	#pragma omp parallel for private(molecular_viscosity)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// calculating the molecular viscosity
		// calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &molecular_viscosity);
		// homogeneous for now
		viscosity_value = min_diff_h_coeff_turb;
		calc_diffusion_coeff(273.15, mean_particle_mass, 1, eff_particle_radius, &molecular_viscosity);
		// the molecular viscosity is the absolute minimum
		if (molecular_viscosity > viscosity_value)
		{
			viscosity_value = molecular_viscosity;
		}
		curl_term_viscosity_eff[i] = viscosity_value;
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











