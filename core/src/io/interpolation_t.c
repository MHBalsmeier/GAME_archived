#include "../enum_and_typedefs.h"
#include "io.h"

int interpolation_t(State *state_0, State *state_p1, State *state_write, double t_0, double t_p1, double t_write)
{
    double weight_0, weight_p1;
    weight_p1 = (t_write - t_0)/(t_p1 - t_0);
    weight_0 = 1 - weight_p1;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_write -> density[i] = weight_0*state_0 -> density[i] + weight_p1*state_p1 -> density[i];
        state_write -> pot_temp[i] = weight_0*state_0 -> pot_temp[i] + weight_p1*state_p1 -> pot_temp[i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        state_write -> wind[i] = weight_0*state_0 -> wind[i] + weight_p1*state_p1 -> wind[i];
    return 1;
}
