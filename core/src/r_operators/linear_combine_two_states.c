#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>
#include <stdio.h>

int linear_combine_two_states(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_out -> density[i] = coeff_0*state_0 -> density[i] + coeff_1*state_1 -> density[i];
        state_out -> density_pot_temp[i] = coeff_0*state_0 -> density_pot_temp[i] + coeff_1*state_1 -> density_pot_temp[i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        state_out -> wind[i] = coeff_0*state_0 -> wind[i] + coeff_1*state_1 -> wind[i];
    return 0;
}
