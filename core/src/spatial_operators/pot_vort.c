/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "../diagnostics/diagnostics.h"

int calc_pot_vort(Vector_field velocity_field, Scalar_field density_field, Curl_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int layer_index, h_index, index, sign, index_for_vertical_gradient, edge_vector_index, edge_vector_index_h, edge_vector_index_dual_area, index_0, index_1, index_2, index_3, sign_0, sign_1, sign_2, sign_3;
    double rhombus_circ, dist_0, dist_1, dist_2, dist_3, dist_0_pre, dist_2_pre, delta_z, covar_0, covar_2, length_rescale_factor, velocity_value, vertical_gradient, density_value;
    for (int i = 0; i < NO_OF_LAYERS*(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H) + NO_OF_DUAL_VECTORS_H; ++i)
    {
        layer_index = i/(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H);
        h_index = i - layer_index*(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H);
        // Rhombus vorticities (stand vertically)
        if (h_index >= NO_OF_DUAL_VECTORS_H)
        {
			edge_vector_index_h = h_index - NO_OF_DUAL_VECTORS_H;
	        edge_vector_index = NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + edge_vector_index_h;
	        edge_vector_index_dual_area = NO_OF_DUAL_VECTORS_H + layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + edge_vector_index_h;
        	rhombus_circ = 0;
        	// The rhombus has four edges.
        	for (int k = 0; k < 4; ++k)
        	{
        		// This is the index of one of the rhombus edges.
				index = NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + k];
			    // This sign is positive for a positive sence of rotation.
			    sign = dualgrid -> vorticity_signs[4*edge_vector_index_h + k];
		    	// This is the value of the velocity at the rhombus edge.
		    	velocity_value = velocity_field[index];
		    	// length_rescale_factor corrects for terrain following coordinates
		    	length_rescale_factor = 1;
		        if (layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
		        {
            		length_rescale_factor = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index]);
            		// In terrain following coordinates, a vertical interpolation onto the reference z coordinate must be made. Therefore, delta_z and vertical_gradient must be determined.
		        	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index];
		        	if (delta_z > 0)
		        		index_for_vertical_gradient = index - NO_OF_VECTORS_PER_LAYER;
		        	else
		        	{
		        		if (layer_index == NO_OF_LAYERS - 1)
		        			index_for_vertical_gradient = index - NO_OF_VECTORS_PER_LAYER;
		        		else
		        			index_for_vertical_gradient = index + NO_OF_VECTORS_PER_LAYER;
		        	}
		        	vertical_gradient = (velocity_field[index] - velocity_field[index_for_vertical_gradient])/(grid -> z_vector[index] - grid -> z_vector[index_for_vertical_gradient]);
		        	// Here, the vertical interpolation is made.
		        	velocity_value += delta_z*vertical_gradient;
    			}
    			rhombus_circ += length_rescale_factor*grid -> normal_distance[index]*sign*velocity_value;
        	}
        	out_field[i] = rhombus_circ/dualgrid -> area[edge_vector_index_dual_area];
        	// We add the result of the curl component of the wind field with the component of the Coriolis vector.
       		out_field[i] += dualgrid -> f_vec[h_index];
       		// The whole stuff must be divided through the density to obtain something similar to a potential vorticity.
			density_value = 0.5*(density_field[layer_index*NO_OF_SCALARS_H + grid -> from_index[edge_vector_index_h]] + density_field[layer_index*NO_OF_SCALARS_H + grid -> to_index[edge_vector_index_h]]);
			out_field[i] = out_field[i]/density_value;
        }
        // tangential vorticities
        else
        {
            if (layer_index == 0 || layer_index == NO_OF_LAYERS)
            {
                index_1 = layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_3 = layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                out_field[i] = 1/dualgrid -> area[layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + h_index]*(grid -> normal_distance[index_1]*sign_1*velocity_field[index_1] + grid -> normal_distance[index_3]*sign_3*velocity_field[index_3]);
       			out_field[i] += dualgrid -> f_vec[h_index];
				density_value = density_field[layer_index*NO_OF_SCALARS_H + h_index];
				if (layer_index == 0)
				{
					density_value = 0.5*(density_field[grid -> from_index[h_index]] + density_field[grid -> to_index[h_index]]);
				}
				else if (layer_index == NO_OF_LAYERS)
				{
					density_value = 0.5*(density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]] + density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]);
				}
				out_field[i] = out_field[i]/density_value;
            }
            else
            {
                index_0 = layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 0] - NO_OF_VECTORS_H;
                index_1 = layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_2 = layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 2] - NO_OF_VECTORS_H;
                index_3 = layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_0 = dualgrid -> h_curl_signs[4*h_index + 0];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_2 = dualgrid -> h_curl_signs[4*h_index + 2];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                dist_0_pre = grid -> normal_distance[index_0];
                dist_1 = grid -> normal_distance[index_1];
                dist_2_pre = grid -> normal_distance[index_2];
                dist_3 = grid -> normal_distance[index_3];
                delta_z = grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> h_curl_indices[4*h_index + 2]]] - grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> h_curl_indices[4*h_index + 2]]];
                dist_0 = sqrt(pow(dist_0_pre, 2) + pow(delta_z, 2));
                delta_z = grid -> z_scalar[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> h_curl_indices[4*h_index + 2]]] - grid -> z_scalar[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> h_curl_indices[4*h_index + 2]]];
                dist_2 = sqrt(pow(dist_2_pre, 2) + pow(delta_z, 2));
                covar_0 = velocity_field[index_0];
                covar_2 = velocity_field[index_2];
                if (layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
                {
                    horizontal_covariant_normalized(velocity_field, layer_index, dualgrid -> h_curl_indices[4*h_index + 2], grid, &covar_0);
                    horizontal_covariant_normalized(velocity_field, layer_index - 1, dualgrid -> h_curl_indices[4*h_index + 2], grid, &covar_2);
                }
                out_field[i] = 1/dualgrid -> area[layer_index*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + h_index]*(dist_0*sign_0*covar_0 + dist_1*sign_1*velocity_field[index_1] + dist_2*sign_2*covar_2 + dist_3*sign_3*velocity_field[index_3]);
       			out_field[i] += dualgrid -> f_vec[h_index];
				density_value = 0.25*(density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]] + density_field[(layer_index - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]] + density_field[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]] + density_field[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]]);
				out_field[i] = out_field[i]/density_value;
            }
        }
    }
    return 0;
}












