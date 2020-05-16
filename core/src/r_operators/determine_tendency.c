#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"
#include "rad.h"
#include "addcomp.h"
#include "surface.h"
#include <stdlib.h>
#include <stdio.h>

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_bool, int add_comps_bool, double delta_t)
{
    Vector_field *density_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density, current_state -> wind, *density_flux, grid);
    Scalar_field *density_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_flux, *density_flux_divergence, grid, 0);
    free(density_flux);
    double mean_particle_mass = M_D/N_A;
    double eff_particle_radius = 130e-12;
    int retval;
    Vector_field *diffusion_mass_flux = malloc(sizeof(Vector_field));
    Scalar_field *mass_diffusion_rate = malloc(sizeof(Scalar_field));
    Scalar_field *temperature = malloc(sizeof(Scalar_field));
    Scalar_field *mass_diffusion_coeff = malloc(sizeof(Scalar_field));
    double mass_diffusion_coeff_para_ratio;
    if (dissipation_on == 1)
    {
        temperature_diagnostics(current_state -> density_pot_temp, current_state -> density, *temperature);
        Vector_field *diffusion_mass_flux_pre = malloc(sizeof(Vector_field));
        retval = grad(current_state -> density, *diffusion_mass_flux_pre, grid);
        for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        {
            calc_diffusion_coeff((*temperature)[i], mean_particle_mass, current_state -> density[i], eff_particle_radius, &(*mass_diffusion_coeff)[i]);
            mass_diffusion_coeff_para_ratio = pow(10, 5);
            (*mass_diffusion_coeff)[i] = mass_diffusion_coeff_para_ratio*(*mass_diffusion_coeff)[i];
        }
        scalar_times_vector(*mass_diffusion_coeff, *diffusion_mass_flux_pre, *diffusion_mass_flux, grid);
        free(diffusion_mass_flux_pre);
        retval = divergence(*diffusion_mass_flux, *mass_diffusion_rate, grid, 0);
        free(diffusion_mass_flux);
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (dissipation_on == 1)
            state_tendency -> density[i] = -(*density_flux_divergence)[i] + (*mass_diffusion_rate)[i] ;
        else
            state_tendency -> density[i] = -(*density_flux_divergence)[i];
    }
    free(density_flux_divergence);
    Vector_field *density_pot_temp_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density_pot_temp, current_state -> wind, *density_pot_temp_flux, grid);
    Scalar_field *density_pot_temp_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_pot_temp_flux, *density_pot_temp_flux_divergence, grid, 0);
    free(density_pot_temp_flux);
    Vector_field *laplace_wind_field = malloc(sizeof(Vector_field));
    laplace_vec(current_state -> wind, *laplace_wind_field, grid, dualgrid);
    Scalar_field *u_dot_laplace_wind = malloc(sizeof(Scalar_field));
    inner(current_state -> wind, *laplace_wind_field, *u_dot_laplace_wind, grid);
    Scalar_field *temp_diffusion_heating = malloc(sizeof(Scalar_field));
    Vector_field *temperature_flux = malloc(sizeof(Vector_field));
    if (dissipation_on == 1)
    {
        Vector_field *temperature_flux_pre = malloc(sizeof(Vector_field));
        retval = grad(*temperature, *temperature_flux_pre, grid);
        scalar_times_vector(*mass_diffusion_coeff, *temperature_flux_pre, *temperature_flux, grid);
        free(temperature_flux_pre);
        free(mass_diffusion_coeff);
        retval = divergence(*temperature_flux, *temp_diffusion_heating, grid, 0);
    }
    free(temperature_flux);
    Scalar_field *pot_temp = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        (*pot_temp)[i] = current_state -> density_pot_temp[i]/current_state -> density[i];
    Scalar_field *exner_pressure = malloc(sizeof(Scalar_field));
    exner_pressure_diagnostics(current_state -> density_pot_temp, *exner_pressure);
    double viscosity_coeff_molecular = 5*pow(10, -5);
    double viscosity_coeff = pow(10, 5)*viscosity_coeff_molecular;
    double friction_heating;
    Scalar_field *rad_heating = calloc(1, sizeof(Scalar_field));
    if (rad_bool == 1)
        retval = calc_rad_heating(*rad_heating, NUMBER_OF_SCALARS);
    double *add_comp_mass_source_rates = calloc(NUMBER_OF_ADD_COMPS*NUMBER_OF_SCALARS, sizeof(double));
    double *add_comp_heat_source_rates = calloc(NUMBER_OF_ADD_COMPS*NUMBER_OF_SCALARS, sizeof(double));
    if (add_comps_bool == 1)
        retval = calc_add_comp_source_rates(add_comp_mass_source_rates, add_comp_heat_source_rates, current_state -> add_comp_densities, current_state -> add_comp_temps, *temperature, NUMBER_OF_ADD_COMPS, NUMBER_OF_SCALARS, delta_t);
    double total_density;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (dissipation_on == 1)
        {
            friction_heating = -viscosity_coeff*current_state -> density[i]*(*u_dot_laplace_wind)[i];
            total_density = current_state -> density[i];
            for (int k = 0; k < NUMBER_OF_ADD_COMPS; ++k)
                total_density += current_state -> add_comp_densities[k*NUMBER_OF_SCALARS + i];
            state_tendency -> density_pot_temp[i] = -(*density_pot_temp_flux_divergence)[i] + (*pot_temp)[i]*(*mass_diffusion_rate)[i] + 1/(C_P*(*exner_pressure)[i])*(friction_heating + (*rad_heating)[i]*current_state -> density[i]/total_density + (*temp_diffusion_heating)[i] + add_comp_heat_source_rates[(NUMBER_OF_ADD_COMPS - 1)*NUMBER_OF_SCALARS + i]);
        }
        else
            state_tendency -> density_pot_temp[i] = -(*density_pot_temp_flux_divergence)[i];
    }
    free(temp_diffusion_heating);
    free(mass_diffusion_rate);
    free(u_dot_laplace_wind);
    free(density_pot_temp_flux_divergence);
    Scalar_field *add_comp_density = malloc(sizeof(Scalar_field));
    Vector_field *add_comp_velocity = malloc(sizeof(Vector_field));
    Vector_field *add_comp_density_flux = malloc(sizeof(Vector_field));
    Scalar_field *add_comp_flux_divergence = malloc(sizeof(Scalar_field));
    int h_index, layer_index;
    for (int i = 0; i < NUMBER_OF_ADD_COMPS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
            (*add_comp_density)[j] = current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j];
        if (i < NUMBER_OF_COND_ADD_COMPS)
        {
            for (int j = 0; j < NUMBER_OF_VECTORS; ++j)
            {
                layer_index = j/NUMBER_OF_VECTORS_PER_LAYER;
                h_index = j - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
                (*add_comp_velocity)[j] = current_state -> wind[j];
                if (h_index < NUMBER_OF_VECTORS_V)
                    (*add_comp_velocity)[j] += ret_sink_velocity(i, 0, 0.001);
            }
            retval = scalar_times_vector(*add_comp_density, *add_comp_velocity, *add_comp_density_flux, grid);
            retval = divergence(*add_comp_density_flux, *add_comp_flux_divergence, grid, 1);
        }
        else
        {
            retval = scalar_times_vector(*add_comp_density, current_state -> wind, *add_comp_density_flux, grid);
            retval = divergence(*add_comp_density_flux, *add_comp_flux_divergence, grid, 0);
        }
        for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
        {
            state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] = -(*add_comp_flux_divergence)[j] + add_comp_mass_source_rates[i*NUMBER_OF_SCALARS + j];
            if (current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j] + delta_t*state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] < 0)
                state_tendency -> add_comp_densities[i*NUMBER_OF_SCALARS + j] = -current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j]/delta_t;
        }
    }
    free(add_comp_mass_source_rates);
    free(add_comp_heat_source_rates);
    free(add_comp_density);
    free(add_comp_density_flux);
    free(add_comp_flux_divergence);
    Scalar_field *add_comp_temp = malloc(sizeof(Scalar_field));
    Scalar_field *add_comp_temp_adv = malloc(sizeof(Scalar_field));
    double c_p_cond;
    for (int i = 0; i < NUMBER_OF_COND_ADD_COMPS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS; ++j)
        {
            layer_index = j/NUMBER_OF_VECTORS_PER_LAYER;
            h_index = j - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
            (*add_comp_velocity)[j] = current_state -> wind[j];
            if (h_index < NUMBER_OF_VECTORS_V)
                (*add_comp_velocity)[j] += ret_sink_velocity(i, 0, 0.001);
        }
        retval = adv_scalar(*add_comp_temp, *add_comp_velocity, *add_comp_temp_adv, grid);
        for (int j = 0; j < NUMBER_OF_SCALARS; ++j)
        {
            c_p_cond = ret_c_p_cond(i, 0, current_state -> add_comp_temps[i*NUMBER_OF_SCALARS + j]);
            total_density = current_state -> density[j];
            for (int k = 0; k < NUMBER_OF_ADD_COMPS; ++k)
                total_density += current_state -> add_comp_densities[k*NUMBER_OF_SCALARS + j];
            if (current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j] > 0)
                state_tendency -> add_comp_temps[i*NUMBER_OF_SCALARS + j] = (*add_comp_temp_adv)[j] + 1/c_p_cond*(*rad_heating)[j]/total_density + 1/(c_p_cond*current_state -> add_comp_densities[i*NUMBER_OF_SCALARS + j])*add_comp_heat_source_rates[i*NUMBER_OF_SCALARS + j];
            else
                state_tendency -> add_comp_temps[i*NUMBER_OF_SCALARS + j] = ((*temperature)[j] - current_state -> add_comp_temps[i*NUMBER_OF_SCALARS + j])/delta_t;
        }
    }
    free(rad_heating);
    free(add_comp_temp);
    free(add_comp_velocity);
    free(add_comp_temp_adv);
    Dual_vector_field *rel_curl = malloc(sizeof(Dual_vector_field));
    curl(current_state -> wind, *rel_curl, grid, dualgrid);
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
    Vector_field *m_pressure_gradient_acc = malloc(sizeof(Vector_field));
    retval = scalar_times_vector(*pot_temp, *exner_pressure_perturb_gradient, *m_pressure_gradient_acc, grid);
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
    int check_energy_bool = 0;
    double enery_integral = 0;
    double z_value, z_g;
    if (check_energy_bool == 1)
    {
        for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        {
            enery_integral += (*e_kin_spec_2)[i]*current_state -> density[i]*grid -> volume[i];
            z_value = grid -> z_scalar[i];
            z_g = z_value/(1 + z_value/RADIUS);
            enery_integral += -grid -> gravity[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_PER_LAYER]*z_g*current_state -> density[i]*grid -> volume[i];
            enery_integral += C_V*(*temperature)[i]*current_state -> density[i]*grid -> volume[i];
        }
        printf("%lf\n", enery_integral);
    }
    free(temperature);
    free(e_kin_spec_2);
    double pot_temp_value, pot_temp_perturb;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (i < NUMBER_OF_VECTORS_V || i >= NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V)
            state_tendency -> wind[i] = 0;
        else
        {
            if (h_index < NUMBER_OF_VECTORS_V)
                pot_temp_value = 0.5*current_state -> density_pot_temp[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H]/current_state -> density[h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H] + 0.5*current_state -> density_pot_temp[h_index + layer_index*NUMBER_OF_SCALARS_H]/current_state -> density[h_index + layer_index*NUMBER_OF_SCALARS_H];
            else
                pot_temp_value = 0.5*(current_state -> density_pot_temp[grid -> to_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H]/current_state -> density[grid -> to_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H] + current_state -> density_pot_temp[grid -> from_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H]/current_state -> density[grid -> from_index[h_index - NUMBER_OF_VECTORS_V] + layer_index*NUMBER_OF_SCALARS_H]);
            pot_temp_perturb = pot_temp_value - grid -> pot_temp_background[i];
            state_tendency -> wind[i] = -(*m_pressure_gradient_acc)[i] + (*abs_curl_tend)[i] - 0.5*(*m_e_kin_tend_2)[i] - pot_temp_perturb/grid -> pot_temp_background[i]*grid -> gravity[i] + dissipation_on*viscosity_coeff*(*laplace_wind_field)[i];
        }
    }
    free(laplace_wind_field);
    free(abs_curl_tend);
    free(m_pressure_gradient_acc);
    free(m_e_kin_tend_2);
    return retval;
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
