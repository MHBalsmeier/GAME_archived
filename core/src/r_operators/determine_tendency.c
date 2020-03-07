#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include "../diagnostics/diagnostics.h"

State tendency(State current_state)
{
    State state_tendency;
    extern Grid grid;
    extern Dualgrid dualgrid;
    Scalar_field adv_density, adv_pot_temp, v_divergence;
    adv_s(current_state.wind, current_state.density, adv_density);
    adv_s(current_state.wind, current_state.pot_temp, adv_pot_temp);
    divergence(current_state.wind, v_divergence);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_tendency.density[i] = -adv_density[i] - current_state.density[i]*v_divergence[i];
        state_tendency.pot_temp[i] = -adv_pot_temp[i];
    }
    Dual_vector_field rel_curl, abs_curl;
    rot(current_state.wind, rel_curl);
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; i++)
        abs_curl[i] = dualgrid.f_vec[i] + rel_curl[i];
    Vector_field abs_curl_tend, pressure_gradient, e_kin_tend_2_m;
    vector_product(current_state.wind, abs_curl, abs_curl_tend);
    Scalar_field pressure, e_kin_spec_2;
    scalar_product(current_state.wind, current_state.wind, e_kin_spec_2);
    grad(e_kin_spec_2, e_kin_tend_2_m);
    pressure_diagnostics(current_state.pot_temp, current_state.density, pressure);
    grad(pressure, pressure_gradient);
    long vert_index, floor_index, h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        vert_index = i/(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        floor_index = vert_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        h_index = i - floor_index;
        if (h_index >= NUMBER_OF_SCALARS_H)
        state_tendency.wind[i] = -(1/current_state.density[i])*pressure_gradient[i] + abs_curl_tend[i] - 0.5*e_kin_tend_2_m[i];
        else
        state_tendency.wind[i] = -(1/current_state.density[i])*pressure_gradient[i] + abs_curl_tend[i] - 0.5*e_kin_tend_2_m[i] + grid.gravity[h_index + vert_index*NUMBER_OF_SCALARS_H];
        if (i < NUMBER_OF_SCALARS_H || i >= NUMBER_OF_VECTORS - NUMBER_OF_SCALARS_H)
            state_tendency.wind[i] = 0;
    }
    return state_tendency;
}
