/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int kinetic_energy(Vector_field in_field, Scalar_field out_field, Grid *grid, int diag_bool)
{
	// It computes only the horizontal kinetic energy. Only this part is neeed for the 3D Lamb transformation.
	// For diag_bool == 3 it computes the 3D kinetic energy. This is only needed for diagnostics.
	int max_index = 6;
	if (diag_bool == 1)
		max_index = 8;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        out_field[i] = 0;
        for (int j = 0; j < max_index; ++j)
        	out_field[i] += grid -> e_kin_weights[8*i + j]*pow(in_field[grid -> e_kin_indices[8*i + j]], 2);
    }
    return 0;
}
