/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include <stdlib.h>
#include <stdio.h>

int linear_combine_two_states(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_out -> density_dry[i] = coeff_0*state_0 -> density_dry[i] + coeff_1*state_1 -> density_dry[i];
        state_out -> entropy_gas[i] = coeff_0*state_0 -> entropy_gas[i] + coeff_1*state_1 -> entropy_gas[i];
        for (int j = 0; j < NUMBER_OF_TRACERS; ++j)
            state_out -> tracer_densities[j*NUMBER_OF_SCALARS + i] = coeff_0*state_0 -> tracer_densities[j*NUMBER_OF_SCALARS + i] + coeff_1*state_1 -> tracer_densities[j*NUMBER_OF_SCALARS + i];
        for (int j = 0; j < NUMBER_OF_CONDENSATED_TRACERS; ++j)
            state_out -> tracer_density_temperatures[j*NUMBER_OF_SCALARS + i] = coeff_0*state_0 -> tracer_density_temperatures[j*NUMBER_OF_SCALARS + i] + coeff_1*state_1 -> tracer_density_temperatures[j*NUMBER_OF_SCALARS + i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        state_out -> velocity_gas[i] = coeff_0*state_0 -> velocity_gas[i] + coeff_1*state_1 -> velocity_gas[i];
    return 0;
}

int linear_combine_two_states_scalars(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        state_out -> density_dry[i] = coeff_0*state_0 -> density_dry[i] + coeff_1*state_1 -> density_dry[i];
        state_out -> entropy_gas[i] = coeff_0*state_0 -> entropy_gas[i] + coeff_1*state_1 -> entropy_gas[i];
        for (int j = 0; j < NUMBER_OF_TRACERS; ++j)
            state_out -> tracer_densities[j*NUMBER_OF_SCALARS + i] = coeff_0*state_0 -> tracer_densities[j*NUMBER_OF_SCALARS + i] + coeff_1*state_1 -> tracer_densities[j*NUMBER_OF_SCALARS + i];
        for (int j = 0; j < NUMBER_OF_CONDENSATED_TRACERS; ++j)
            state_out -> tracer_density_temperatures[j*NUMBER_OF_SCALARS + i] = coeff_0*state_0 -> tracer_density_temperatures[j*NUMBER_OF_SCALARS + i] + coeff_1*state_1 -> tracer_density_temperatures[j*NUMBER_OF_SCALARS + i];
    }
    return 0;
}
