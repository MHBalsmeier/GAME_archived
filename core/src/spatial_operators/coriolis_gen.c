/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int coriolis_gen(Vector_field a_field, Curl_field b_field, Vector_field out_field, Grid *grid)
{
    double component_0, component_1, component_2, component_3, term_0, term_1;
    int layer_index, h_index, retval;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            retval = trsk_modified(a_field, b_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &term_0, grid);
            retval = recov_hor_par_curl(b_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &component_1, grid);
            retval = recov_hor_ver_pri(a_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &component_2, grid);
        }
        else
        {
            retval = recov_ver_0_pri(a_field, layer_index, h_index, &component_0, grid);
            retval = recov_ver_0_curl(b_field, layer_index, h_index, &component_1, grid);
            retval = recov_ver_1_pri(a_field, layer_index, h_index, &component_2, grid);
            retval = recov_ver_1_curl(b_field, layer_index, h_index, &component_3, grid);
        	term_0 = component_0*component_3;
        }
        term_1 = component_1*component_2;
        out_field[i] = term_0 - term_1;
    }
    return retval;
}
