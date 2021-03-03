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
    double upper_weight, lower_weight, upper_value, lower_value, z_lower, z_upper;
	#pragma omp parallel for private(layer_index, h_index, i, upper_weight, lower_weight, upper_value, lower_value, z_lower, z_upper)
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
        	Horizontal non-traditional component (horizontal potential vorticity times vertical mass flux density).
            -------------------------------------------------------------------------------------------------------
            */
            // z position of the upper end of the face
            z_upper = 0.5*(
            grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index - NO_OF_SCALARS_H]]
            + grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index - NO_OF_SCALARS_H]]);
            // z position of the lower end of the face
            z_lower = 0.5*(
            grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index - NO_OF_SCALARS_H]]
            + grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index - NO_OF_SCALARS_H]]);
            // vorticity flux value at upper end of the face
            upper_value = 0.5*(
            mass_flux_density[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index - NO_OF_SCALARS_H]]
            + mass_flux_density[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index - NO_OF_SCALARS_H]])
            *pot_vorticity[h_index - NO_OF_SCALARS_H + layer_index*2*NO_OF_VECTORS_H];
            // vorticity flux value at lower end of the face
            lower_value = 0.5*(
            mass_flux_density[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index - NO_OF_SCALARS_H]]
            + mass_flux_density[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index - NO_OF_SCALARS_H]])
            *pot_vorticity[h_index - NO_OF_SCALARS_H  + (layer_index + 1)*2*NO_OF_VECTORS_H];
            // determining the weights
            upper_weight = find_volume(pow((RADIUS + grid -> z_vector[i])/(RADIUS + z_lower), 2), RADIUS + grid -> z_vector[i], RADIUS + z_upper)/find_volume(1, RADIUS + z_lower, RADIUS + z_upper);
            lower_weight = find_volume(1, RADIUS + z_lower, RADIUS + grid -> z_vector[i])/find_volume(1, RADIUS + z_lower, RADIUS + z_upper);
            // adding to the result
            out_field[i] += -(upper_weight*upper_value + lower_weight*lower_value);
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









