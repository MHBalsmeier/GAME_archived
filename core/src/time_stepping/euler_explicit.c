#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int euler_explicit(State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid)
{
    State *tendency_0 = malloc(sizeof(State));
    tendency(state_0, tendency_0, grid, dualgrid);
    linear_combine_two_states(state_0, tendency_0, state_p1, 1, delta_t);
    free(tendency_0);
    return 0;
}
