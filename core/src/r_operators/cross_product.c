#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>
#include <stdio.h>

int cross_product(Vector_field a_field, Dual_vector_field b_field, Vector_field out_field, Grid *grid)
{
    double component_0, component_1, component_2, component_3, term_0, term_1;
    long layer_index, h_index;
    int retval;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - NUMBER_OF_VECTORS_PER_LAYER*layer_index;
        if(h_index >= NUMBER_OF_VECTORS_V)
        {
            retval = recov_hor_par_pri(a_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &component_0, grid);
            retval = recov_hor_par_dual(b_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &component_1, grid);
            retval = recov_hor_ver_pri(a_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &component_2, grid);
            retval = recov_hor_ver_dual(b_field, layer_index, h_index - NUMBER_OF_VECTORS_V, &component_3, grid);
        }
        else
        {
            retval = recov_ver_0_pri(a_field, layer_index, h_index, &component_0, grid);
            retval = recov_ver_0_dual(b_field, layer_index, h_index, &component_1, grid);
            retval = recov_ver_1_pri(a_field, layer_index, h_index, &component_2, grid);
            retval = recov_ver_1_dual(b_field, layer_index, h_index, &component_3, grid);
        }
        term_0 = component_0*component_3;
        term_1 = component_2*component_1;
        out_field[i] = term_0 - term_1;
    }
    return 0;
}
