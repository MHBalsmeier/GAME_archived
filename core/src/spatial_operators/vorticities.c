/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/
/*
Here, vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "geos95.h"
#include "../thermodynamics/thermodynamics.h"
#include "spatial_operators.h"

int calc_pot_vort(Vector_field velocity_field, Scalar_field density_field, Diagnostics *diagnostics, Grid *grid, Dualgrid *dualgrid)
{
	// It is called "potential vorticity", but it is not Ertel's potential vorticity. It is the absolute vorticity divided by the density.
	calc_rel_vort(velocity_field, diagnostics -> rel_vort, grid, dualgrid);
	// pot_vort is a misuse of name here
	add_f_to_rel_vort(diagnostics -> rel_vort, diagnostics -> pot_vort, dualgrid);
    int i, layer_index, h_index, edge_vector_index_h, upper_from_index, upper_to_index;
    double density_value;
    // determining the density value by which we need to divide
    #pragma omp parallel for private (i, layer_index, h_index, edge_vector_index_h, upper_from_index, upper_to_index, density_value)
    for (i = 0; i < NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H; ++i)
    {
        layer_index = i/(2*NO_OF_VECTORS_H);
        h_index = i - layer_index*2*NO_OF_VECTORS_H;
        
        /*
        Determination of the density value by which we have to divide.
        */
        // interpolation of the density to the center or the rhombi
        if (h_index >= NO_OF_VECTORS_H)
        {
			edge_vector_index_h = h_index - NO_OF_VECTORS_H;
			density_value = 0;
			for (int j = 0; j < 4; ++j)
			{
				density_value
				+= grid -> density_to_rhombus_weights[4*edge_vector_index_h + j]
				*density_field[layer_index*NO_OF_SCALARS_H + grid -> density_to_rhombus_indices[4*edge_vector_index_h + j]];
			}
        }
        // interpolation of the density to the half level edges
        else
        {
        	// linear extrapolation to the TOA
        	if (layer_index == 0)
        	{
				density_value
				= 0.5*(density_field[grid -> from_index[h_index]] + density_field[grid -> to_index[h_index]])
				// the gradient
				+ (0.5*(density_field[grid -> from_index[h_index]] + density_field[grid -> to_index[h_index]])
				- 0.5*(density_field[grid -> from_index[h_index] + NO_OF_SCALARS_H] + density_field[grid -> to_index[h_index] + NO_OF_SCALARS_H]))
				/(grid -> z_vector[NO_OF_SCALARS + h_index] - grid -> z_vector[NO_OF_SCALARS + NO_OF_VECTORS_PER_LAYER + h_index])
				// delta z
				*(grid -> z_vector[0] - grid -> z_vector[NO_OF_SCALARS + h_index]);
        	}
        	// linear extrapolation to the surface
            else if (layer_index == NO_OF_LAYERS)
            {
				density_value =
				0.5*(density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]] + density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]])
				// the gradient
				+ (0.5*(density_field[(layer_index - 2)*NO_OF_SCALARS_H + grid -> from_index[h_index]] + density_field[(layer_index - 2)*NO_OF_SCALARS_H + grid -> to_index[h_index]])
				- 0.5*(density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]] + density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]))
				/(grid -> z_vector[NO_OF_SCALARS + (layer_index - 2)*NO_OF_VECTORS_PER_LAYER + h_index] - grid -> z_vector[NO_OF_SCALARS + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + h_index])
				// delta z
				*(0.5*(grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]] + grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]])
				- grid -> z_vector[NO_OF_SCALARS + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + h_index]);
            }
            else
            {
            	upper_from_index = (layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index];
            	upper_to_index = (layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index];
            	density_value = 0.25*(density_field[upper_from_index] + density_field[upper_to_index]
            	+ density_field[upper_from_index + NO_OF_SCALARS_H] + density_field[upper_to_index + NO_OF_SCALARS_H]);
            }
        }
        
        // division by the density to obtain the "potential vorticity"
		diagnostics -> pot_vort[i] = diagnostics -> pot_vort[i]/density_value;
    }
    return 0;
}

int add_f_to_rel_vort(Curl_field rel_vort, Curl_field out_field, Dualgrid *dualgrid)
{
	/*
	Adding the Coriolis parameter to the relative vorticity.
	*/
    int i, layer_index, h_index;
    #pragma omp parallel for private(i, layer_index, h_index)
    for (i = 0; i < NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H; ++i)
    {
        layer_index = i/(2*NO_OF_VECTORS_H);
        h_index = i - layer_index*2*NO_OF_VECTORS_H;
   		out_field[i] = rel_vort[i] + dualgrid -> f_vec[h_index];
    }
    return 0;
}

int calc_rel_vort(Vector_field velocity_field, Curl_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int i, layer_index, h_index, index, sign, index_for_vertical_gradient, edge_vector_index, edge_vector_index_h, index_0, index_1, index_2, index_3;
    double rhombus_circ, dist_0, dist_1, dist_2, dist_3, delta_z, covar_0, covar_2, length_rescale_factor, velocity_value, vertical_gradient;
	#pragma omp parallel for private(i, layer_index, h_index, index, sign, index_for_vertical_gradient, edge_vector_index, edge_vector_index_h, index_0, index_1, index_2, index_3, rhombus_circ, dist_0, dist_1, dist_2, dist_3, delta_z, covar_0, covar_2, length_rescale_factor, velocity_value, vertical_gradient)
    for (i = NO_OF_VECTORS_H; i < NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H; ++i)
    {
        layer_index = i/(2*NO_OF_VECTORS_H);
        h_index = i - layer_index*2*NO_OF_VECTORS_H;
        // Rhombus vorticities (stand vertically)
        if (h_index >= NO_OF_VECTORS_H)
        {
			edge_vector_index_h = h_index - NO_OF_VECTORS_H;
	        edge_vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + edge_vector_index_h;
        	rhombus_circ = 0;
        	// The rhombus has four edges.
        	for (int k = 0; k < 4; ++k)
        	{
        		// This is the index of one of the rhombus edges.
				index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + k];
			    // This sign is positive for a positive sence of rotation.
			    sign = dualgrid -> vorticity_signs[4*edge_vector_index_h + k];
		    	// This is the value of the velocity at the rhombus edge.
		    	velocity_value = velocity_field[index];
		    	// length_rescale_factor corrects for terrain following coordinates
		    	length_rescale_factor = 1;
		        if (layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
		        {
            		length_rescale_factor = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index]);
            		// In terrain following coordinates, a vertical interpolation onto the reference z coordinate must be made. Therefore, delta_z and vertical_gradient must be determined.
		        	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index];
		        	if (delta_z > 0)
		        	{
		        		index_for_vertical_gradient = index - NO_OF_VECTORS_PER_LAYER;
	        		}
		        	else
		        	{
		        		if (layer_index == NO_OF_LAYERS - 1)
		        		{
		        			index_for_vertical_gradient = index - NO_OF_VECTORS_PER_LAYER;
	        			}
		        		else
		        		{
		        			index_for_vertical_gradient = index + NO_OF_VECTORS_PER_LAYER;
	        			}
		        	}
		        	vertical_gradient = (velocity_field[index] - velocity_field[index_for_vertical_gradient])/(grid -> z_vector[index] - grid -> z_vector[index_for_vertical_gradient]);
		        	// Here, the vertical interpolation is made.
		        	velocity_value += delta_z*vertical_gradient;
    			}
    			rhombus_circ += length_rescale_factor*grid -> normal_distance[index]*sign*velocity_value;
        	}
        	out_field[i] = rhombus_circ/dualgrid -> area[i];
        }
        // tangential vorticities
        else
        {
        	// At the lower boundary, w vanishes. Furthermore, the covariant velocity at the surface is also zero.
            if (layer_index == NO_OF_LAYERS)
            {
                index_2 = layer_index*NO_OF_VECTORS_PER_LAYER - NO_OF_VECTORS_H + h_index;
                dist_2 = grid -> normal_distance[index_2];
                horizontal_covariant(velocity_field, layer_index - 1, h_index, grid, &covar_2);
                out_field[i] = 1/dualgrid -> area[layer_index*2*NO_OF_VECTORS_H + h_index]*dist_2*covar_2;
            }
            else
            {
                index_0 = layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_SCALARS_H + h_index;
                index_1 = layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index];
                index_2 = layer_index*NO_OF_VECTORS_PER_LAYER - NO_OF_VECTORS_H + h_index;
                index_3 = layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index];
                dist_0 = grid -> normal_distance[index_0];
                dist_1 = grid -> normal_distance[index_1];
                dist_2 = grid -> normal_distance[index_2];
                dist_3 = grid -> normal_distance[index_3];
                horizontal_covariant(velocity_field, layer_index, h_index, grid, &covar_0);
                horizontal_covariant(velocity_field, layer_index - 1, h_index, grid, &covar_2);
                out_field[i] = 1/dualgrid -> area[i]*(-dist_0*covar_0 + dist_1*velocity_field[index_1] + dist_2*covar_2 - dist_3*velocity_field[index_3]);
            }
        }
    }
    // At the upper boundary, the tangential vorticity is assumed to have no vertical shear.
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	out_field[i] = out_field[i + 2*NO_OF_VECTORS_H];
    }
    return 0;
}









