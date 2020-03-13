#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid)
{
    Scalar_field *adv_density = malloc(sizeof(Scalar_field));
    Scalar_field *adv_density_pot_temp = malloc(sizeof(Scalar_field));
    Scalar_field *v_divergence = malloc(sizeof(Scalar_field));
    adv_s(current_state -> wind, current_state -> density, *adv_density, grid);
    adv_s(current_state -> wind, current_state -> density_pot_temp, *adv_density_pot_temp, grid);
    divergence(current_state -> wind, *v_divergence, grid);
    Vector_field *density_pot_temp_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density_pot_temp, current_state -> wind, *density_pot_temp_flux, grid);
    Scalar_field *density_pot_temp_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_pot_temp_flux, *density_pot_temp_flux_divergence, grid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_tendency -> density[i] = (*adv_density)[i] - (current_state -> density[i])*(*v_divergence)[i];
        state_tendency -> density_pot_temp[i] = (*density_pot_temp_flux_divergence)[i];
    }
    free(adv_density);
    free(adv_density_pot_temp);
    free(density_pot_temp_flux_divergence);
    free(v_divergence);
    free(density_pot_temp_flux);
    Dual_vector_field *rel_curl = malloc(sizeof(Dual_vector_field));
    Dual_vector_field *abs_curl = malloc(sizeof(Dual_vector_field));
    rot(current_state -> wind, *rel_curl, grid, dualgrid);
    long layer_index;
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; i++)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        (*abs_curl)[i] = dualgrid -> f_vec[i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_VECTORS_V)] + (*rel_curl)[i];
    }
    free(rel_curl);
    Vector_field *abs_curl_tend = malloc(sizeof(Vector_field));
    Vector_field *e_kin_tend_2_m = malloc(sizeof(Vector_field));
    vector_product(current_state -> wind, *abs_curl, *abs_curl_tend, grid);
    free(abs_curl);
    Scalar_field *exner_pressure = malloc(sizeof(Scalar_field));
    Scalar_field *e_kin_spec_2 = malloc(sizeof(Scalar_field));
    scalar_product(current_state -> wind, current_state -> wind, *e_kin_spec_2, grid);
    grad(*e_kin_spec_2, *e_kin_tend_2_m, grid);
    free(e_kin_spec_2);
    exner_pressure_diagnostics(current_state -> density_pot_temp, *exner_pressure);
    Vector_field *exner_pressure_gradient = malloc(sizeof(Vector_field));
    grad(*exner_pressure, *exner_pressure_gradient, grid);
    free(exner_pressure);
    Scalar_field *pot_temp = malloc(sizeof(Scalar_field));
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        (*pot_temp)[i] = current_state -> density_pot_temp[i]/current_state -> density[i];
    Vector_field *m_pressure_gradient_acc = malloc(sizeof(Vector_field));
    scalar_times_vector(*pot_temp, *exner_pressure_gradient, *m_pressure_gradient_acc, grid);
    free(pot_temp);
    free(exner_pressure_gradient);
    double density;
    long h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (i < NUMBER_OF_VECTORS_V || i >= NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V)
            state_tendency -> wind[i] = 0;
        else
            state_tendency -> wind[i] = -(*m_pressure_gradient_acc)[i] + (*abs_curl_tend)[i] - 0.5*(*e_kin_tend_2_m)[i] + grid -> gravity[i];
    }
    free(abs_curl_tend);
    free(m_pressure_gradient_acc);
    free(e_kin_tend_2_m);
    return 0;
}
