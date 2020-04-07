#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid)
{
    Vector_field *density_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density, current_state -> wind, *density_flux, grid);
    Scalar_field *density_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_flux, *density_flux_divergence, grid);
    free(density_flux);
    double mass_molecular_diffusion_coeff, mass_diffusion_coeff;
    double mean_particle_mass = M_D/N_A;
    double eff_particle_radius = 130e-12;
    Scalar_field *laplace_density = malloc(sizeof(Scalar_field));
    laplace(current_state -> density, *laplace_density, grid);
    Scalar_field *exner_pressure = malloc(sizeof(Scalar_field));
    exner_pressure_diagnostics(current_state -> density_pot_temp, *exner_pressure);
    Scalar_field *temperature = malloc(sizeof(Scalar_field));
    temperature_diagnostics(current_state -> density_pot_temp, current_state -> density, *temperature);
    double mean_temperature = 250;
    double mean_density = 0.7;
    int retval = calc_diffusion_coeff(mean_temperature, mean_particle_mass, mean_density, eff_particle_radius, &mass_molecular_diffusion_coeff);
    mass_diffusion_coeff = pow(10, 5)*mass_molecular_diffusion_coeff;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        state_tendency -> density[i] = -(*density_flux_divergence)[i] + mass_diffusion_coeff*(*laplace_density)[i];
    free(density_flux_divergence);
    Vector_field *density_pot_temp_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density_pot_temp, current_state -> wind, *density_pot_temp_flux, grid);
    Scalar_field *density_pot_temp_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_pot_temp_flux, *density_pot_temp_flux_divergence, grid);
    free(density_pot_temp_flux);
    Vector_field *laplace_wind_field = malloc(sizeof(Vector_field));
    laplace_vec(current_state -> wind, *laplace_wind_field, grid, dualgrid);
    Scalar_field *u_dot_friction = malloc(sizeof(Scalar_field));
    inner(current_state -> wind, *laplace_wind_field, *u_dot_friction, grid);
    double heat_power_density;
    Scalar_field *laplace_temperature = malloc(sizeof(Scalar_field));
    laplace(*temperature, *laplace_temperature, grid);
    double temperature_molecular_diffusion_coeff, temperature_diffusion_coeff;
    retval = calc_diffusion_coeff(mean_temperature, mean_particle_mass, mean_density, eff_particle_radius, &temperature_molecular_diffusion_coeff);
    temperature_diffusion_coeff = pow(10, 5)*temperature_molecular_diffusion_coeff;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        heat_power_density = -current_state -> density[i]*(*u_dot_friction)[i];
        heat_power_density += temperature_diffusion_coeff*(*laplace_temperature)[i];
        state_tendency -> density_pot_temp[i] = -(*density_pot_temp_flux_divergence)[i] + current_state -> density_pot_temp[i]/current_state -> density[i]*mass_diffusion_coeff*(*laplace_density)[i] + 1/(C_P*(*exner_pressure)[i])*heat_power_density;
    }
    free(laplace_temperature);
    free(temperature);
    free(u_dot_friction);
    free(laplace_density);
    free(density_pot_temp_flux_divergence);
    Dual_vector_field *rel_curl = malloc(sizeof(Dual_vector_field));
    curl(current_state -> wind, *rel_curl, grid, dualgrid);
    long layer_index;
    Dual_vector_field *abs_curl = malloc(sizeof(Dual_vector_field));
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        (*abs_curl)[i] = dualgrid -> f_vec[i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER] + (*rel_curl)[i];
    }
    free(rel_curl);
    Scalar_field *exner_pressure_perturb = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        (*exner_pressure_perturb)[i] = (*exner_pressure)[i] - grid -> exner_pressure_background[i];
    free(exner_pressure);
    Vector_field *exner_pressure_perturb_gradient = malloc(sizeof(Vector_field));
    grad(*exner_pressure_perturb, *exner_pressure_perturb_gradient, grid);
    free(exner_pressure_perturb);
    Scalar_field *pot_temp = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        (*pot_temp)[i] = current_state -> density_pot_temp[i]/current_state -> density[i];
    Vector_field *m_pressure_gradient_acc = malloc(sizeof(Vector_field));
    scalar_times_vector(*pot_temp, *exner_pressure_perturb_gradient, *m_pressure_gradient_acc, grid);
    free(pot_temp);
    free(exner_pressure_perturb_gradient);
    Vector_field *abs_curl_tend = malloc(sizeof(Vector_field));
    cross_product(current_state -> wind, *abs_curl, *abs_curl_tend, grid);
    free(abs_curl);
    Scalar_field *e_kin_spec_2 = malloc(sizeof(Scalar_field));
    inner(current_state -> wind, current_state -> wind, *e_kin_spec_2, grid);
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        (*m_pressure_gradient_acc)[i] = C_P*(*m_pressure_gradient_acc)[i];
    Vector_field *m_e_kin_tend_2 = malloc(sizeof(Vector_field));
    grad(*e_kin_spec_2, *m_e_kin_tend_2, grid);
    free(e_kin_spec_2);
    double viscosity_coeff = 5;
    double pot_temp_perturb;
    long h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (i < NUMBER_OF_VECTORS_V || i >= NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V)
            state_tendency -> wind[i] = 0;
        else
        {
            state_tendency -> wind[i] = -(*m_pressure_gradient_acc)[i] + (*abs_curl_tend)[i] - 0.5*(*m_e_kin_tend_2)[i] + viscosity_coeff*(*laplace_wind_field)[i];
            if (h_index < NUMBER_OF_VECTORS_V)
            {
                pot_temp_perturb = 0.5*current_state -> density_pot_temp[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H]/current_state -> density[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H] + 0.5*current_state -> density_pot_temp[h_index + layer_index*NUMBER_OF_SCALARS_H]/current_state -> density[h_index + layer_index*NUMBER_OF_SCALARS_H] - grid -> pot_temp_background[h_index + layer_index*NUMBER_OF_VECTORS_V];
                state_tendency -> wind[i] -= C_P*pot_temp_perturb*grid -> exner_pressure_background_gradient[layer_index*NUMBER_OF_VECTORS_V + h_index];
            }
        }
    }
    free(laplace_wind_field);
    free(abs_curl_tend);
    free(m_pressure_gradient_acc);
    free(m_e_kin_tend_2);
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
