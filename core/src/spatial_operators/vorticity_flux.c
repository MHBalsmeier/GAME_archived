/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int vorticity_flux(Vector_field mass_flux_density, Curl_field pot_vorticity, Vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int layer_index, h_index, i;
	#pragma omp parallel for private(layer_index, h_index, i)
    for (i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        
        /*
        Calculating the horizontal component of the vorticity flux term.
        ----------------------------------------------------------------
        */
        if (h_index >= NO_OF_SCALARS_H)
        {
        	/*
        	Traditional component (vertical potential vorticity times horizontal mass flux density).
            ----------------------------------------------------------------------------------------
            */
            vorticity_flux_horizontal_traditional(mass_flux_density, pot_vorticity, layer_index, h_index - NO_OF_SCALARS_H, &out_field[i], grid);
            
        	/*
        	Horizontal "non-standard" component (horizontal potential vorticity times vertical mass flux density).
            -------------------------------------------------------------------------------------------------------
            */
            out_field[i]
			= out_field[i]
			- 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index - NO_OF_SCALARS_H]) + 6]
			*mass_flux_density[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index - NO_OF_SCALARS_H]]
			*pot_vorticity[h_index - NO_OF_SCALARS_H + layer_index*2*NO_OF_VECTORS_H];
			out_field[i]
			= out_field[i]
			- 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index - NO_OF_SCALARS_H]) + 7]
			*mass_flux_density[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index - NO_OF_SCALARS_H]]
			*pot_vorticity[h_index - NO_OF_SCALARS_H + (layer_index + 1)*2*NO_OF_VECTORS_H];
			out_field[i]
			= out_field[i]
			- 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index - NO_OF_SCALARS_H]) + 6]
			*mass_flux_density[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index - NO_OF_SCALARS_H]]
			*pot_vorticity[h_index - NO_OF_SCALARS_H + layer_index*2*NO_OF_VECTORS_H];
			out_field[i]
			= out_field[i]
			- 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index - NO_OF_SCALARS_H]) + 7]
			*mass_flux_density[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index - NO_OF_SCALARS_H]]
			*pot_vorticity[h_index - NO_OF_SCALARS_H + (layer_index + 1)*2*NO_OF_VECTORS_H];
        }
        
        /*
        Calculating the vertical component of the vorticity flux term.
        --------------------------------------------------------------
        */
        else
        {    
			vorticity_flux_vertical(mass_flux_density, pot_vorticity, layer_index, h_index, &out_field[i], grid, dualgrid);
        }
    }
    return 0;
}









