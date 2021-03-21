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

int calc_horizontal_shear(State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function calculates the shear of the horizontal wind field at vector points.
	*/
	// calculating eastward and northward winds at edge
	calc_uv_at_edge(state -> velocity_gas, diagnostics -> u_at_edge, diagnostics -> v_at_edge, grid);
	// averaging the wind components to the cell centers
	edges_to_cells(diagnostics -> u_at_edge, diagnostics -> u_at_cell, grid);
	edges_to_cells(diagnostics -> v_at_edge, diagnostics -> v_at_cell, grid);
	// computing the horizontal gradient of the eastward wind
	// u_at_edge and v_at_edge are misusages of name from now on
	grad_hor(diagnostics -> u_at_cell, diagnostics -> u_at_edge, grid);
	// computing the horizontal gradient of the northward wind
	grad_hor(diagnostics -> v_at_cell, diagnostics -> v_at_edge, grid);
	int layer_index, h_index;
	double comp_orth, comp_tang, dudx, dudy, dvdx, dvdy;
	#pragma omp parallel for private(layer_index, h_index, dudx, dudy, dvdx, dvdy, comp_orth, comp_tang)
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		// only the shear of the horizontal wind field is diagnozed
		if (h_index >= NO_OF_SCALARS_H)
		{
			// diagnozing u quantities
			comp_orth = diagnostics -> u_at_edge[layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
			tangential_wind(diagnostics -> u_at_edge, layer_index, h_index - NO_OF_SCALARS_H, &comp_tang, grid);
			passive_turn(comp_orth, comp_tang, -grid -> direction[h_index - NO_OF_SCALARS_H], &dudx, &dudy);
			// diagnozing v quantities
			comp_orth = diagnostics -> v_at_edge[layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
			tangential_wind(diagnostics -> v_at_edge, layer_index, h_index - NO_OF_SCALARS_H, &comp_tang, grid);
			passive_turn(comp_orth, comp_tang, -grid -> direction[h_index - NO_OF_SCALARS_H], &dvdx, &dvdy);
			// calculating the deformation according to the MPAS paper
			diagnostics -> shear[i] = sqrt(pow(dudx - dvdy, 2) + pow(dudy + dvdx, 2));
		}
	}
	return 0;
}
