/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int coriolis_gen(Vector_field a_field, Curl_field b_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index;
    double component_0, component_1, component_2, component_3, term_0, term_1;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
        	/*
        	Only one term is needed here, do not doubt. See Lamb transformation in spherical coordinates (only the horizontal velocity is used in the cross product).
            */
            trsk_modified(a_field, b_field, layer_index, h_index - NO_OF_SCALARS_H, &term_0, grid);
		    out_field[i] = term_0;
        }
        else
        {
        	// Here, two terms are indeed required.
            recov_ver_0_pri(a_field, layer_index, h_index, &component_0, grid);
            recov_ver_0_curl(b_field, layer_index, h_index, &component_1, grid);
            recov_ver_1_pri(a_field, layer_index, h_index, &component_2, grid);
            recov_ver_1_curl(b_field, layer_index, h_index, &component_3, grid);
        	term_0 = component_0*component_3;
		    term_1 = component_1*component_2;
		    out_field[i] = term_0 - term_1;
        }
    }
    return 0;
}
