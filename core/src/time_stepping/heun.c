#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int heun(State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid)
{
    State *tendency_0 = malloc(sizeof(State));
    tendency(state_0, tendency_0, grid, dualgrid);
    State *state_star = malloc(sizeof(State));
    linear_combine_two_states(state_0, tendency_0, state_star, 1, delta_t);
    State *tendency_star = malloc(sizeof(State));
    tendency(state_star, tendency_star, grid, dualgrid);
    free(state_star);
    State *tendency_eff = malloc(sizeof(State));
    linear_combine_two_states(tendency_0, tendency_star, tendency_eff, 0.5, 0.5);
    free(tendency_star);
    free(tendency_0);
    linear_combine_two_states(state_0, tendency_eff, state_p1, 1, delta_t);
    free(tendency_eff);
    return 0;
}
