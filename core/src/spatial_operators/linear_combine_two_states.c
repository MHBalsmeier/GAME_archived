/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include <stdlib.h>
#include <stdio.h>

int linear_combine_two_states(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1)
{
	int i, j;
	#pragma omp parallel for private(i, j)
    for (i = 0; i < NO_OF_SCALARS; ++i)
    {
        state_out -> temperature_gas[i] = coeff_0*state_0 -> temperature_gas[i] + coeff_1*state_1 -> temperature_gas[i];
        
        for (j = 0; j < NO_OF_CONSTITUENTS; ++j)
        {
            state_out -> mass_densities[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> mass_densities[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> mass_densities[j*NO_OF_SCALARS + i];
        }
        
        for (j = 0; j < NO_OF_CONDENSED_CONSTITUENTS; ++j)
        {
            state_out -> condensed_density_temperatures[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> condensed_density_temperatures[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> condensed_density_temperatures[j*NO_OF_SCALARS + i];
        }
        
        for (j = 0; j < NO_OF_CONSTITUENTS; ++j)
        {
            state_out -> entropy_densities[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> entropy_densities[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> entropy_densities[j*NO_OF_SCALARS + i];
        }
    }
    
    #pragma omp parallel for private(i)
    for (i = 0; i < NO_OF_VECTORS; ++i)
    {
        state_out -> velocity_gas[i] = coeff_0*state_0 -> velocity_gas[i] + coeff_1*state_1 -> velocity_gas[i];
    }
    return 0;
}










