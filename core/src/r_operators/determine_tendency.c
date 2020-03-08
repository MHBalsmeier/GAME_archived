#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int tendency(State *current_state, State *state_tendency, Grid *grid, Dualgrid *dualgrid)
{
    Scalar_field *adv_density = malloc(sizeof(Scalar_field));
    Scalar_field *adv_pot_temp = malloc(sizeof(Scalar_field));
    Scalar_field *v_divergence = malloc(sizeof(Scalar_field));
    adv_s(current_state -> wind, current_state -> density, *adv_density, grid);
    adv_s(current_state -> wind, current_state -> pot_temp, *adv_pot_temp, grid);
    divergence(current_state -> wind, *v_divergence, grid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_tendency -> density[i] = -*adv_density[i] - current_state -> density[i]**v_divergence[i];
        state_tendency -> pot_temp[i] = -*adv_pot_temp[i];
    }
    free(adv_density);
    free(adv_pot_temp);
    free(v_divergence);
    Dual_vector_field *rel_curl = malloc(sizeof(Dual_vector_field));
    Dual_vector_field *abs_curl = malloc(sizeof(Dual_vector_field));
    rot(current_state -> wind, *rel_curl, grid, dualgrid);
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; i++)
        (*abs_curl)[0] = dualgrid -> f_vec[i] + (*rel_curl)[0];
    free(rel_curl);
    Vector_field *abs_curl_tend = malloc(sizeof(Vector_field));
    Vector_field *pressure_gradient = malloc(sizeof(Vector_field));
    Vector_field *e_kin_tend_2_m = malloc(sizeof(Vector_field));
    vector_product(current_state -> wind, *abs_curl, *abs_curl_tend, grid);
    free(abs_curl);
    Scalar_field *pressure = malloc(sizeof(Scalar_field));
    Scalar_field *e_kin_spec_2 = malloc(sizeof(Scalar_field));
    scalar_product(current_state -> wind, current_state -> wind, *e_kin_spec_2, grid);
    grad(*e_kin_spec_2, *e_kin_tend_2_m, grid);
    free(e_kin_spec_2);
    pressure_diagnostics(current_state -> pot_temp, current_state -> density, *pressure);
    grad(*pressure, *pressure_gradient, grid);
    free(pressure);
    long vert_index, floor_index, h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        vert_index = i/(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        floor_index = vert_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        h_index = i - floor_index;
        if (h_index >= NUMBER_OF_SCALARS_H)
        state_tendency -> wind[i] = -(1/current_state -> density[i])*(*pressure_gradient)[i] + (*abs_curl_tend)[i] - 0.5*(*e_kin_tend_2_m)[i];
        else
        state_tendency -> wind[i] = -(1/current_state -> density[i])*(*pressure_gradient)[i] + (*abs_curl_tend)[i] - 0.5*(*e_kin_tend_2_m)[i] + grid -> gravity[h_index + vert_index*NUMBER_OF_SCALARS_H];
        if (i < NUMBER_OF_SCALARS_H || i >= NUMBER_OF_VECTORS - NUMBER_OF_SCALARS_H)
            state_tendency -> wind[i] = 0;
    }
    free(abs_curl_tend);
    free(pressure_gradient);
    free(e_kin_tend_2_m);
    return 0;
}
