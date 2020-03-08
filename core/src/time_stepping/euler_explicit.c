#include "../enum_and_typedefs.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int euler_explicit(State *state_0, State *tendency_0, State *state_p1, double delta_t)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_p1 -> density[i] = state_0 -> density[i] + delta_t*tendency_0 -> density[i];
        state_p1 -> pot_temp[i] = state_0 -> pot_temp[i] + delta_t*tendency_0 -> pot_temp[i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        state_p1 -> wind[i] = state_0 -> wind[i] + delta_t*tendency_0 -> wind[i];
    Scalar_field *pressure = malloc(sizeof(Scalar_field));
    pressure_diagnostics(state_p1 -> pot_temp, state_p1 -> density, *pressure);
    free(pressure);
    return 1;
}
