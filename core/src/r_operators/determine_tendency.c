#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid)
{
    Vector_field *density_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density, current_state -> wind, *density_flux, grid);
    Scalar_field *density_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_flux, *density_flux_divergence, grid);
    free(density_flux);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        state_tendency -> density[i] = -(*density_flux_divergence)[i];
    free(density_flux_divergence);
    Vector_field *density_pot_temp_flux = malloc(sizeof(Vector_field));
    scalar_times_vector(current_state -> density_pot_temp, current_state -> wind, *density_pot_temp_flux, grid);
    Scalar_field *density_pot_temp_flux_divergence = malloc(sizeof(Scalar_field));
    divergence(*density_pot_temp_flux, *density_pot_temp_flux_divergence, grid);
    free(density_pot_temp_flux);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        state_tendency -> density_pot_temp[i] = -(*density_pot_temp_flux_divergence)[i];
    free(density_pot_temp_flux_divergence);
    Dual_vector_field *rel_curl = malloc(sizeof(Dual_vector_field));
    curl(current_state -> wind, *rel_curl, grid, dualgrid);
    long layer_index;
    Dual_vector_field *abs_curl = malloc(sizeof(Dual_vector_field));
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; i++)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        (*abs_curl)[i] = dualgrid -> f_vec[i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER] + (*rel_curl)[i];
    }
    free(rel_curl);
    free(abs_curl);
    Scalar_field *e_kin_spec_2 = malloc(sizeof(Scalar_field));
    scalar_product(current_state -> wind, current_state -> wind, *e_kin_spec_2, grid);
    Vector_field *m_e_kin_tend_2 = malloc(sizeof(Vector_field));
    grad(*e_kin_spec_2, *m_e_kin_tend_2, grid);
    free(e_kin_spec_2);
    Scalar_field *exner_pressure = malloc(sizeof(Scalar_field));
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
    Vector_field *abs_curl_tend = malloc(sizeof(Vector_field));
    vector_product(current_state -> wind, *abs_curl, *abs_curl_tend, grid);
    double density;
    long h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        (*m_pressure_gradient_acc)[i] = C_P*(*m_pressure_gradient_acc)[i];
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (i < NUMBER_OF_VECTORS_V || i >= NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V)
            state_tendency -> wind[i] = 0;
        else
            state_tendency -> wind[i] = -(*m_pressure_gradient_acc)[i] + (*abs_curl_tend)[i] - 0.5*(*m_e_kin_tend_2)[i] + grid -> gravity[i];
    }
    free(abs_curl_tend);
    free(m_pressure_gradient_acc);
    free(m_e_kin_tend_2);
    return 0;
}
