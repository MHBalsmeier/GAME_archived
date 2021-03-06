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

int hori_viscosity_eff(State *state, Scalar_field viscosity_eff, Grid *grid, Diagnostics *diagnostics, Forcings *forcings, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	double molecular_viscosity;
	// the maximum diffusion coefficient (based on stability, including a 30 % safety margin)
	double max_diff_h_coeff_turb = (1 - 0.3)*0.25*grid -> mean_area_edge/delta_t;
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_edge*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	double viscosity_value, shear;
	inner_product(forcings -> rel_vort_tend, forcings -> rel_vort_tend, diagnostics -> scalar_field_placeholder, grid, 0);
	#pragma omp parallel for private(molecular_viscosity, shear)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// this is only an estimate
		shear = 0;
		if (diagnostics -> e_kin[i] > EPSILON_SECURITY)
		{
			shear = sqrt(diagnostics -> scalar_field_placeholder[i])/sqrt(diagnostics -> e_kin[i]);
		}
		// preliminary result
		viscosity_value =  grid -> mean_area_edge*config_info -> diff_h_smag_fac*shear;
		
		/*
		Checking if the calculated value is "allowed".
		*/
		// calculating the molecular viscosity
		calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &molecular_viscosity);
		// the molecular viscosity is the absolute minimum
		if (viscosity_value < molecular_viscosity)
		{
			viscosity_value = molecular_viscosity;
		}
		// turbulent minimum
		if (viscosity_value < min_diff_h_coeff_turb)
		{
			viscosity_value = min_diff_h_coeff_turb;
		}
		// stability
		if (viscosity_value > max_diff_h_coeff_turb)
		{
			viscosity_value = max_diff_h_coeff_turb;
		}
		
		// result
		viscosity_eff[i] = viscosity_value;
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











