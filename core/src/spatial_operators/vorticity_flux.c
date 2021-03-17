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
    int layer_index, h_index, h_index_shifted, i, number_of_edges, j;
    double vert_weight;
	#pragma omp parallel for private(layer_index, h_index, i, number_of_edges, vert_weight, j, h_index_shifted)
    for (i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        out_field[i] = 0;
        
        /*
        Calculating the horizontal component of the vorticity flux term.
        ----------------------------------------------------------------
        */
        if (h_index >= NO_OF_SCALARS_H)
        {
        	h_index_shifted = h_index - NO_OF_SCALARS_H;
        	/*
        	Traditional component (vertical potential vorticity times horizontal mass flux density).
            ----------------------------------------------------------------------------------------
            */
			// From_index comes before to_index as usual.
			if (grid -> from_index[h_index_shifted] < NO_OF_PENTAGONS)
			{
				for (j = 0; j < 4; ++j)
				{
					out_field[i] +=
					grid -> trsk_weights[10*h_index_shifted + j]
					*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index_shifted + j]]
					*pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
				}
			}
			else
			{
				for (j = 0; j < 5; ++j)
				{
					if (j == 2)
					{
						out_field[i] +=
						grid -> trsk_weights[10*h_index_shifted + j]
						*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index_shifted + j]]
						*0.5
						*(pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]]
						+ pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + h_index_shifted]);
					}
					else
					{
						out_field[i] +=
						grid -> trsk_weights[10*h_index_shifted + j]
						*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index_shifted + j]]
						*pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
					}
				}
			}
			if (grid -> to_index[h_index_shifted] < NO_OF_PENTAGONS)	
			{
				for (j = 5; j < 9; ++j)
				{
					out_field[i] +=
					grid -> trsk_weights[10*h_index_shifted + j]
					*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index_shifted + j]]
					*pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
				}
			}
			else
			{
				for (j = 5; j < 10; ++j)
				{
					if (j == 7)
					{
						out_field[i] +=
						grid -> trsk_weights[10*h_index_shifted + j]
						*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index_shifted + j]]
						*0.5
						*(pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]]
						+ pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + h_index_shifted]);
					}
					else
					{
						out_field[i] +=
						grid -> trsk_weights[10*h_index_shifted + j]
						*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index_shifted + j]]
						*pot_vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
					}
				}
			}
            
        	/*
        	Horizontal "non-standard" component (horizontal potential vorticity times vertical mass flux density).
            -------------------------------------------------------------------------------------------------------
            */
            out_field[i]
			-= 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index_shifted]) + 6]
			*mass_flux_density[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index_shifted]]
			*pot_vorticity[h_index_shifted + layer_index*2*NO_OF_VECTORS_H];
			out_field[i]
			-= 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index_shifted]) + 7]
			*mass_flux_density[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index_shifted]]
			*pot_vorticity[h_index_shifted + (layer_index + 1)*2*NO_OF_VECTORS_H];
			out_field[i]
			-= 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index_shifted]) + 6]
			*mass_flux_density[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index_shifted]]
			*pot_vorticity[h_index_shifted + layer_index*2*NO_OF_VECTORS_H];
			out_field[i]
			-= 0.5
			*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index_shifted]) + 7]
			*mass_flux_density[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index_shifted]]
			*pot_vorticity[h_index_shifted + (layer_index + 1)*2*NO_OF_VECTORS_H];
        }
        
        /*
        Calculating the vertical component of the vorticity flux term.
        --------------------------------------------------------------
        */
        else
        {    
			/*
			Determining the vertical acceleration due to the vorticity flux term.
			*/
			// determining the number of edges
			number_of_edges = 6;
			if (h_index < NO_OF_PENTAGONS)
			{
				number_of_edges = 5;
			}
			// determining the vertical interpolation weight
			vert_weight = 0.5;
			if (layer_index == 0 || layer_index == NO_OF_LAYERS)
			{
				vert_weight = 1;
			}
			if (layer_index >= 1)
			{
				for (j = 0; j < number_of_edges; ++j)
				{
					out_field[i] +=
					vert_weight
					*grid -> inner_product_weights[8*((layer_index - 1)*NO_OF_SCALARS_H + h_index) + j]
					*mass_flux_density[NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
					*pot_vorticity[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
				}
			}
			if (layer_index <= NO_OF_LAYERS - 1)
			{
				for (j = 0; j < number_of_edges; ++j)
				{
					out_field[i] +=
					vert_weight
					*grid -> inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + h_index) + j]
					*mass_flux_density[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
					*pot_vorticity[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
				}
			}
        }
    }
    return 0;
}









