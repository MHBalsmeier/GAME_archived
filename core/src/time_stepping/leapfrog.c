#include "../enum_and_typedefs.h"
#include "time_stepping.h"
#include "../diagnostics/diagnostics.h"

State leapfrog(State state_m1, State state_0, State tendency_0, double delta_t)
{
    State state_p1;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_p1.density[i] = state_m1.density[i] + 2*delta_t*tendency_0.density[i];
        state_p1.pot_temp[i] = state_m1.pot_temp[i] + 2*delta_t*tendency_0.pot_temp[i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        state_p1.wind[i] = state_m1.wind[i] + 2*delta_t*tendency_0.wind[i];
    Scalar_field pressure;
    pressure_diagnostics(state_p1.pot_temp, state_p1.density, pressure);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        state_p1.pressure[i] = pressure[i];
    return state_p1;
}
