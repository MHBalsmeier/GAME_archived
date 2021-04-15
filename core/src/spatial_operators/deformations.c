/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include "../spatial_operators/spatial_operators.h"
#include "geos95.h"
#include <stdlib.h>
#include <stdio.h>

int calc_deformations(State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function calculates the deformations needed for the horizontal momentum diffusion.
	*/
	
	// diagonal component
	divv_h(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> deform_diag[i] = 5.0/3*diagnostics -> velocity_gas_divv[i];
	}
	
	// off-diagonal component
	int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		diagnostics -> deform_off_diag[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index]
		= -diagnostics -> rel_vort[NO_OF_VECTORS_H + 2*layer_index*NO_OF_VECTORS_H + h_index];
	}
	return 0;
}
