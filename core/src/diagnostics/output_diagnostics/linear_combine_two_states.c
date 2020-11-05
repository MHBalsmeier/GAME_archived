/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include <stdlib.h>
#include <stdio.h>

int linear_combine_two_states(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1)
{
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        state_out -> density_dry[i] = coeff_0*state_0 -> density_dry[i] + coeff_1*state_1 -> density_dry[i];
        state_out -> entropy_density_dry[i] = coeff_0*state_0 -> entropy_density_dry[i] + coeff_1*state_1 -> entropy_density_dry[i];
        state_out -> temperature_gas[i] = coeff_0*state_0 -> temperature_gas[i] + coeff_1*state_1 -> temperature_gas[i];
        for (int j = 0; j < NO_OF_TRACERS; ++j)
        {
            state_out -> tracer_densities[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> tracer_densities[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> tracer_densities[j*NO_OF_SCALARS + i];
        }
        for (int j = 0; j < NO_OF_CONDENSED_TRACERS; ++j)
        {
            state_out -> tracer_density_temperatures[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> tracer_density_temperatures[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> tracer_density_temperatures[j*NO_OF_SCALARS + i];
        }
        for (int j = 0; j < NO_OF_GASEOUS_TRACERS; ++j)
        {
            state_out -> tracer_entropy_densities[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> tracer_entropy_densities[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> tracer_entropy_densities[j*NO_OF_SCALARS + i];
        }
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        state_out -> velocity_gas[i] = coeff_0*state_0 -> velocity_gas[i] + coeff_1*state_1 -> velocity_gas[i];
    }
    return 0;
}









