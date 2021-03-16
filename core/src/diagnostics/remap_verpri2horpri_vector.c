/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This function remaps vertical primal vectors to horizontal primal vector points.
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int remap_verpri2horpri_vector(Vector_field vector_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component
    = grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + 6]
    *vector_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    *component
    += grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + 7]
    *vector_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]];
    *component
    += grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]) + 6]
    *vector_field[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
    *component
    += grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]) + 7]
    *vector_field[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]];
    *component = 0.5*(*component);
    return 0;
}
