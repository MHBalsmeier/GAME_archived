#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int runge_kutta_third_order(State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid, int dissipation_on, int rad_on, int add_comps_on)
{
    State *state_star = malloc(sizeof(State));
    State *tendency_0 = malloc(sizeof(State));
    tendency(state_0, tendency_0, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    linear_combine_two_states(state_0, tendency_0, state_star, 1, delta_t/3.0);
    State *state_sub_1_mod = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_0, state_sub_1_mod, 1, delta_t);
    free(tendency_0);
    linear_combine_two_states(state_sub_1_mod, state_0, state_star, 1.0/3.0, 2.0/3.0);
    free(state_sub_1_mod);
    State *tendency_star = malloc(sizeof(State));
    tendency(state_star, tendency_star, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    free(state_star);
    State *state_star_star = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_star, state_star_star, 1, delta_t/2.0);
    State *state_sub_2_mod = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_star, state_sub_2_mod, 1, delta_t);
    free(tendency_star);
    State *tendency_star_star = malloc(sizeof(State));
    tendency(state_star_star, tendency_star_star, grid, dualgrid, dissipation_on, rad_on, add_comps_on, delta_t);
    linear_combine_two_states(state_sub_2_mod, state_0, state_star_star, 0.5, 0.5);
    free(state_sub_2_mod);
    free(state_star_star);
    linear_combine_two_states(state_0, tendency_star_star, state_p1, 1, delta_t);
    free(tendency_star_star);
    return 0;
}
