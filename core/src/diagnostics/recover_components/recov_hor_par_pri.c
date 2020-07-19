/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
This is a recovery function at horizontal vector points. If layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS, a vertical interpolation of the input field onto the desired z coordinate must be done.
*/

#include "../../enum_and_typedefs.h"

int recov_hor_par_pri(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
	for (int i = 0; i < 10; ++i)
		*component += grid -> trsk_modified_weights[10*h_index + i]*in_field[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_modified_velocity_indices[10*h_index + i]];
    return 0;
}
