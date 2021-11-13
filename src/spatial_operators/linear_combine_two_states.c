/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains a function for linearly combining two states.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"

int linear_combine_two_states(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1, Grid *grid)
{
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        
        for (int j = 0; j < NO_OF_CONSTITUENTS; ++j)
        {
            state_out -> rho[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> rho[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> rho[j*NO_OF_SCALARS + i];
        }
        
        state_out -> rhotheta[i] = coeff_0*state_0 -> rhotheta[i] + coeff_1*state_1 -> rhotheta[i];
        state_out -> theta_pert[i] = state_out -> rhotheta[i]/state_out -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] - grid -> theta_bg[i];
        state_out -> exner_pert[i] = coeff_0*state_0 -> exner_pert[i] + coeff_1*state_1 -> exner_pert[i];
        
        for (int j = 0; j < NO_OF_CONDENSED_CONSTITUENTS; ++j)
        {
            state_out -> condensed_density_temperatures[j*NO_OF_SCALARS + i]
            = coeff_0*state_0 -> condensed_density_temperatures[j*NO_OF_SCALARS + i]
            + coeff_1*state_1 -> condensed_density_temperatures[j*NO_OF_SCALARS + i];
        }
    }
    
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        state_out -> wind[i] = coeff_0*state_0 -> wind[i] + coeff_1*state_1 -> wind[i];
    }
    return 0;
}










