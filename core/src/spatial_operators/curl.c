/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "../diagnostics/diagnostics.h"

int curl(Vector_field in_field, Curl_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int layer_index, h_index, index, sign, retval, index_for_vertical_gradient, edge_vector_index, edge_vector_index_h, edge_vector_index_dual_area;
    double rhombus_circ, dist_0, dist_1, dist_2, dist_3, dist_0_pre, dist_2_pre, delta_z, covar_0, covar_2, length_rescale_factor, velocity_value, vertical_gradient;
    int index_0, index_1, index_2, index_3, sign_0, sign_1, sign_2, sign_3;
    for (int i = 0; i < NUMBER_OF_LAYERS*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H) + NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        layer_index = i/(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
        h_index = i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
        {
			edge_vector_index_h = h_index - NUMBER_OF_DUAL_VECTORS_H;
	        edge_vector_index = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + edge_vector_index_h;
	        edge_vector_index_dual_area = NUMBER_OF_DUAL_VECTORS_H + layer_index*(NUMBER_OF_VECTORS_H + NUMBER_OF_DUAL_VECTORS_H) + edge_vector_index_h;
        	rhombus_circ = 0;
        	for (int k = 0; k < 4; ++k)
        	{
				index = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + k];
			    sign = dualgrid -> vorticity_signs[4*edge_vector_index_h + k];
			    if (sign != -1 && sign != 1)
			    	printf("Error in curl, position 0.\n");
		    	velocity_value = in_field[index];
		    	length_rescale_factor = 1;
		        if (layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
		        {
            		length_rescale_factor = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index]);
		        	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index];
		        	if (delta_z > 0)
		        		index_for_vertical_gradient = index - NUMBER_OF_VECTORS_PER_LAYER;
		        	else
		        	{
		        		if (layer_index == NUMBER_OF_LAYERS - 1)
		        			index_for_vertical_gradient = index - NUMBER_OF_VECTORS_PER_LAYER;
		        		else
		        			index_for_vertical_gradient = index + NUMBER_OF_VECTORS_PER_LAYER;
		        	}
		        	vertical_gradient = (in_field[index] - in_field[index_for_vertical_gradient])/(grid -> z_vector[index] - grid -> z_vector[index_for_vertical_gradient]);
		        	velocity_value += delta_z*vertical_gradient;
    			}
    			rhombus_circ += length_rescale_factor*grid -> normal_distance[index]*sign*velocity_value;
        	}
        	out_field[i] = rhombus_circ/dualgrid -> area[edge_vector_index_dual_area];
        }
        else
        {
            if (layer_index == 0 || layer_index == NUMBER_OF_LAYERS)
            {
                index_1 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_3 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                out_field[i] = 1/dualgrid -> area[layer_index*(NUMBER_OF_VECTORS_H + NUMBER_OF_DUAL_VECTORS_H) + h_index]*(grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_3]*sign_3*in_field[index_3]);
            }
            else
            {
                index_0 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 0] - NUMBER_OF_VECTORS_H;
                index_1 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_2 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 2] - NUMBER_OF_VECTORS_H;
                index_3 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_0 = dualgrid -> h_curl_signs[4*h_index + 0];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_2 = dualgrid -> h_curl_signs[4*h_index + 2];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                dist_0_pre = grid -> normal_distance[index_0];
                dist_1 = grid -> normal_distance[index_1];
                dist_2_pre = grid -> normal_distance[index_2];
                dist_3 = grid -> normal_distance[index_3];
                delta_z = grid -> z_scalar[layer_index*NUMBER_OF_SCALARS_H + grid -> to_index[dualgrid -> h_curl_indices[4*h_index + 2]]] - grid -> z_scalar[layer_index*NUMBER_OF_SCALARS_H + grid -> from_index[dualgrid -> h_curl_indices[4*h_index + 2]]];
                dist_0 = sqrt(pow(dist_0_pre, 2) + pow(delta_z, 2));
                delta_z = grid -> z_scalar[(layer_index - 1)*NUMBER_OF_SCALARS_H + grid -> to_index[dualgrid -> h_curl_indices[4*h_index + 2]]] - grid -> z_scalar[(layer_index - 1)*NUMBER_OF_SCALARS_H + grid -> from_index[dualgrid -> h_curl_indices[4*h_index + 2]]];
                dist_2 = sqrt(pow(dist_2_pre, 2) + pow(delta_z, 2));
                covar_0 = in_field[index_0];
                covar_2 = in_field[index_2];
                if (layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
                {
                    retval = horizontal_covariant_normalized(in_field, layer_index, dualgrid -> h_curl_indices[4*h_index + 2], grid, &covar_0);
                    if (retval != 0)
                        printf("Error in horizontal_covariant_normalized called at position 0 from curl.\n");
                    retval = horizontal_covariant_normalized(in_field, layer_index - 1, dualgrid -> h_curl_indices[4*h_index + 2], grid, &covar_2);
                    if (retval != 0)
                        printf("Error in horizontal_covariant_normalized called at position 1 from curl.\n");
                }
                out_field[i] = 1/dualgrid -> area[layer_index*(NUMBER_OF_VECTORS_H + NUMBER_OF_DUAL_VECTORS_H) + h_index]*(dist_0*sign_0*covar_0 + dist_1*sign_1*in_field[index_1] + dist_2*sign_2*covar_2 + dist_3*sign_3*in_field[index_3]);
            }
        }
    }
    return 0;
}












