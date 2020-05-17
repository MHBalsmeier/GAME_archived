#include "../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int calc_temp_diffusion_coeffs(Scalar_field *temperature, Scalar_field *density, Scalar_field *mass_diffusion_coeff_numerical_h, Scalar_field *mass_diffusion_coeff_numerical_v)
{
    double mean_particle_mass = M_D/N_A;
    double eff_particle_radius = 130e-12;
    double mass_diffusion_coeff, mass_diffusion_coeff_para_ratio_h, mass_diffusion_coeff_para_ratio_v;
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        calc_diffusion_coeff((*temperature)[i], mean_particle_mass, (*density)[i], eff_particle_radius, &mass_diffusion_coeff);
        mass_diffusion_coeff_para_ratio_h = pow(10, 5);
		mass_diffusion_coeff_para_ratio_v = pow(10, 4);
        (*mass_diffusion_coeff_numerical_h)[i] = mass_diffusion_coeff_para_ratio_h*mass_diffusion_coeff;
        (*mass_diffusion_coeff_numerical_v)[i] = mass_diffusion_coeff_para_ratio_v*mass_diffusion_coeff;
    }
	return 0;
}

int calc_mass_diffusion_coeffs(Scalar_field *temperature, Scalar_field *density, Scalar_field *temp_diffusion_coeff_numerical_h, Scalar_field *temp_diffusion_coeff_numerical_v)
{
    double mean_particle_mass = M_D/N_A;
    double eff_particle_radius = 130e-12;
    double temp_diffusion_coeff, temp_diffusion_coeff_para_ratio_h, temp_diffusion_coeff_para_ratio_v;
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        calc_diffusion_coeff((*temperature)[i], mean_particle_mass, (*density)[i], eff_particle_radius, &temp_diffusion_coeff);
        temp_diffusion_coeff_para_ratio_h = pow(10, 5);
		temp_diffusion_coeff_para_ratio_v = pow(10, 4);
        (*temp_diffusion_coeff_numerical_h)[i] = temp_diffusion_coeff_para_ratio_h*temp_diffusion_coeff;
        (*temp_diffusion_coeff_numerical_v)[i] = temp_diffusion_coeff_para_ratio_v*temp_diffusion_coeff;
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
