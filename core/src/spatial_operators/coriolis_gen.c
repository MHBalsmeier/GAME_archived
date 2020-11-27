/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int vorticity_flux(Vector_field a_field, Curl_field b_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            vorticity_flux_horizontal(a_field, b_field, layer_index, h_index - NO_OF_SCALARS_H, &out_field[i], grid);
        }
        else
        {    
			vorticity_flux_vertical(a_field, b_field, layer_index, h_index, &out_field[i], grid);
        }
    }
    return 0;
}
