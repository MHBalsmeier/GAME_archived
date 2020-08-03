/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>
#include <omp.h>

int kinetic_energy(Vector_field in_field, Scalar_field out_field, Grid *grid, int diag_bool)
{
	// It computes only the horizontal kinetic energy. Only this part is neeed for the 3D Lamb transformation.
	// For diag_bool == 1 it computes the 3D kinetic energy. This is only needed for diagnostics.
	int layer_index, h_index;
	int i, j;
	#pragma omp parallel for private (i, j, layer_index, h_index)
    for (i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
        out_field[i] = 0;
        for (j = 0; j < 6; ++j)
        {
        	out_field[i] += grid -> e_kin_weights[8*i + j]*pow(in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]], 2);
    	}
    	if (diag_bool == 1)
    	{
        	out_field[i] += grid -> e_kin_weights[8*i + 6]*pow(in_field[h_index + layer_index*NO_OF_VECTORS_PER_LAYER], 2);
        	out_field[i] += grid -> e_kin_weights[8*i + 7]*pow(in_field[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER], 2);
    	}
    }
    return 0;
}
