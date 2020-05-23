#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int runge_kutta_fourth_order(State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_on, int add_comps_on)
{
    State *tendency_0 = malloc(sizeof(State));
    tendency(state_0, tendency_0, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    State *state_half_from_below = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_0, state_half_from_below, 1, delta_t/2.0);
    State *tendency_half_from_below = malloc(sizeof(State));
    tendency(state_half_from_below, tendency_half_from_below, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    free(state_half_from_below);
    State *state_half_from_above = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_half_from_below, state_half_from_above, 1, delta_t/2.0);
    State *tendency_half_from_above = malloc(sizeof(State));
    tendency(state_half_from_above, tendency_half_from_above, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    free(state_half_from_above);
    State *state_star = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_half_from_below, state_star, 1, delta_t);
    State *tendency_star = malloc(sizeof(State));
    tendency(state_star, tendency_star, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    free(state_star);
    State *tendency_eff_a = malloc(sizeof(State));
    linear_combine_two_states(tendency_0, tendency_half_from_below, tendency_eff_a, 1.0/6.0, 2.0/6.0);
    free(tendency_half_from_below);
    free(tendency_0);
    State *tendency_eff_b = malloc(sizeof(State));
    linear_combine_two_states(tendency_half_from_above, tendency_star, tendency_eff_b, 2.0/6.0, 1.0/6.0);
    free(tendency_half_from_above);
    free(tendency_star);
    State *tendency_eff = malloc(sizeof(State));
    linear_combine_two_states(tendency_eff_a, tendency_eff_b, tendency_eff, 1, 1);
    free(tendency_eff_a);
    free(tendency_eff_b);
    linear_combine_two_states(state_0, tendency_eff, state_p1, 1, delta_t);
    free(tendency_eff);
    return 0;
}
