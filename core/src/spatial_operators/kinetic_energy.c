/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int kinetic_energy(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	// It computes the horizontal kinetic energy.
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        out_field[i] = 0;
        for (int j = 0; j < 6; ++j)
        	out_field[i] += grid -> e_kin_weights[6*i + j]*pow(in_field[grid -> e_kin_indices[6*i + j]], 2);
    }
    return 0;
}
