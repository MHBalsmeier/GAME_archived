#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int adams_bashforth(State *state_m1, State *state_0, State *state_p1, double delta_t, Grid *grid, Dualgrid *dualgrid)
{
    State *tendency_m1 = malloc(sizeof(State));
    State *tendency_0 = malloc(sizeof(State));
    tendency(state_m1, tendency_m1, grid, dualgrid);
    tendency(state_0, tendency_0, grid, dualgrid);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_p1 -> density[i] = state_m1 -> density[i] + delta_t*(1.5*tendency_0 -> density[i] - 0.5*tendency_m1 -> density[i]);
        state_p1 -> density_pot_temp[i] = state_m1 -> density_pot_temp[i] + delta_t*(1.5*tendency_0 -> density_pot_temp[i] - 0.5*tendency_m1 -> density_pot_temp[i]);
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        state_p1 -> wind[i] = state_m1 -> wind[i] + delta_t*(1.5*tendency_0 -> wind[i] - 0.5*tendency_m1 -> wind[i]);
    free(tendency_m1);
    free(tendency_0);
    return 0;
}
