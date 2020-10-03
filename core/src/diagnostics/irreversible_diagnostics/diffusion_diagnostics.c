/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>


int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int calc_mass_diffusion_coeffs(State *state, Config_info *config_info, Scalar_field mass_diffusion_coeff_numerical_h, Scalar_field mass_diffusion_coeff_numerical_v)
{
	if (config_info -> mass_diff_h == 1 || config_info -> mass_diff_v == 1)
	{
		double mean_particle_mass = M_D/N_A;
		double eff_particle_radius = 130e-12;
		double mass_diffusion_coeff, mass_diffusion_coeff_para_ratio_h, mass_diffusion_coeff_para_ratio_v;
		#pragma omp parallel for private (mass_diffusion_coeff, mass_diffusion_coeff_para_ratio_h, mass_diffusion_coeff_para_ratio_v)
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> density_dry[i], eff_particle_radius, &mass_diffusion_coeff);
		    if (config_info -> mass_diff_h == 1)
		    {
		    	mass_diffusion_coeff_para_ratio_h = pow(10, 5);
	    	}
		    else
		    {
		    	mass_diffusion_coeff_para_ratio_h = 0;
	    	}
		    if (config_info -> mass_diff_v == 1)
		    {
		    	mass_diffusion_coeff_para_ratio_v = pow(10, 3);
	    	}
		    else
		    {
		    	mass_diffusion_coeff_para_ratio_v = 0;
	    	}
		    mass_diffusion_coeff_numerical_h[i] = mass_diffusion_coeff_para_ratio_h*mass_diffusion_coeff;
		    mass_diffusion_coeff_numerical_v[i] = mass_diffusion_coeff_para_ratio_v*mass_diffusion_coeff;
		}
	}
	return 0;
}

int calc_temp_diffusion_coeffs(State *state, Config_info *config_info, Scalar_field temp_diffusion_coeff_numerical_h, Scalar_field temp_diffusion_coeff_numerical_v)
{
	if (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1)
	{
		double mean_particle_mass = M_D/N_A;
		double eff_particle_radius = 130e-12;
		double temp_diffusion_coeff, temp_diffusion_coeff_para_ratio_h, temp_diffusion_coeff_para_ratio_v, rho_h, c_h_v;
		#pragma omp parallel for private (temp_diffusion_coeff, temp_diffusion_coeff_para_ratio_h, temp_diffusion_coeff_para_ratio_v)
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> density_dry[i], eff_particle_radius, &temp_diffusion_coeff);
		    if (config_info -> temperature_diff_h == 1)
		    {
				temp_diffusion_coeff_para_ratio_h = pow(10, 5);
	    	}
			else
		    {
				temp_diffusion_coeff_para_ratio_h = 0;
	    	}
		    if (config_info -> temperature_diff_v == 1)
		    {
				temp_diffusion_coeff_para_ratio_v = pow(10, 3);
	    	}
			else
		    {
				temp_diffusion_coeff_para_ratio_v = 0;
	    	}
			rho_h = state -> density_dry[i] + state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
			c_h_v = spec_heat_cap_diagnostics_v(state -> density_dry[i], state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
		    temp_diffusion_coeff_numerical_h[i] = temp_diffusion_coeff_para_ratio_h*rho_h*c_h_v*temp_diffusion_coeff;
		    temp_diffusion_coeff_numerical_v[i] = temp_diffusion_coeff_para_ratio_v*rho_h*c_h_v*temp_diffusion_coeff;
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
