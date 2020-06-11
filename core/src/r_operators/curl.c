#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "../diagnostics/diagnostics.h"

int curl(Vector_field in_field, Dual_vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int layer_index, h_index, index_0, index_1, index_2, index_3, sign_0, sign_1, sign_2, sign_3, retval, index_for_vertical_gradient, edge_vector_index, edge_vector_index_h, edge_vector_index_dual_area;
    double rhombus_circ, dist_0, dist_1, dist_2, dist_3, dist_0_pre, dist_2_pre, delta_z, covar_0, covar_2, length_rescale_factor_0, length_rescale_factor_1, length_rescale_factor_2, length_rescale_factor_3, velocity_value_0, velocity_value_1, velocity_value_2, velocity_value_3, vertical_gradient;
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
        {
            out_field[i] = 0;
            for (int j = 0; j < 3; ++j)
            {
            	edge_vector_index_h = dualgrid -> adjacent_vector_indices_h[3*(h_index - NUMBER_OF_DUAL_VECTORS_H) + j];
                index_0 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + 0];
                index_1 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + 1];
                index_2 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + 2];
                index_3 = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices[4*edge_vector_index_h + 3];
                sign_0 = dualgrid -> vorticity_signs[4*edge_vector_index_h + 0];
                if (sign_0 != -1 && sign_0 != 1)
                	printf("Error in curl, position 0.\n");
                sign_1 = dualgrid -> vorticity_signs[4*edge_vector_index_h + 1];
                if (sign_1 != -1 && sign_1 != 1)
                	printf("Error in curl, position 1.\n");
                sign_2 = dualgrid -> vorticity_signs[4*edge_vector_index_h + 2];
                if (sign_2 != -1 && sign_2 != 1)
                	printf("Error in curl, position 2.\n");
                sign_3 = dualgrid -> vorticity_signs[4*edge_vector_index_h + 3];
                if (sign_3 != -1 && sign_3 != 1)
                	printf("Error in curl, position 3.\n");
                edge_vector_index = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + edge_vector_index_h;
                edge_vector_index_dual_area = NUMBER_OF_DUAL_VECTORS_H + layer_index*(NUMBER_OF_VECTORS_H + NUMBER_OF_DUAL_VECTORS_H) + edge_vector_index_h;
                if (layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
                {
                	length_rescale_factor_0 = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index_0]);
                	length_rescale_factor_1 = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index_1]);
                	length_rescale_factor_2 = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index_2]);
                	length_rescale_factor_3 = (RADIUS + grid -> z_vector[edge_vector_index])/(RADIUS + grid -> z_vector[index_3]);
                	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index_0];
                	if (delta_z > 0)
                		index_for_vertical_gradient = index_0 - NUMBER_OF_VECTORS_PER_LAYER;
                	else
                	{
                		if (layer_index == NUMBER_OF_LAYERS - 1)
                			index_for_vertical_gradient = index_0 - NUMBER_OF_VECTORS_PER_LAYER;
                		else
                			index_for_vertical_gradient = index_0 + NUMBER_OF_VECTORS_PER_LAYER;
                	}
                	vertical_gradient = (in_field[index_0] - in_field[index_for_vertical_gradient])/(grid -> z_vector[index_0] - grid -> z_vector[index_for_vertical_gradient]);
                	velocity_value_0 = in_field[index_0] + delta_z*vertical_gradient;
                	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index_1];
                	if (delta_z > 0)
                		index_for_vertical_gradient = index_1 - NUMBER_OF_VECTORS_PER_LAYER;
                	else
                	{
                		if (layer_index == NUMBER_OF_LAYERS - 1)
                			index_for_vertical_gradient = index_1 - NUMBER_OF_VECTORS_PER_LAYER;
                		else
                			index_for_vertical_gradient = index_1 + NUMBER_OF_VECTORS_PER_LAYER;
                	}
                	vertical_gradient = (in_field[index_1] - in_field[index_for_vertical_gradient])/(grid -> z_vector[index_1] - grid -> z_vector[index_for_vertical_gradient]);
                	velocity_value_1 = in_field[index_1] + delta_z*vertical_gradient;
                	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index_2];
                	if (delta_z > 0)
                		index_for_vertical_gradient = index_2 - NUMBER_OF_VECTORS_PER_LAYER;
                	else
                	{
                		if (layer_index == NUMBER_OF_LAYERS - 1)
                			index_for_vertical_gradient = index_2 - NUMBER_OF_VECTORS_PER_LAYER;
                		else
                			index_for_vertical_gradient = index_2 + NUMBER_OF_VECTORS_PER_LAYER;
                	}
                	vertical_gradient = (in_field[index_2] - in_field[index_for_vertical_gradient])/(grid -> z_vector[index_2] - grid -> z_vector[index_for_vertical_gradient]);
                	velocity_value_2 = in_field[index_2] + delta_z*vertical_gradient;
                	delta_z = grid -> z_vector[edge_vector_index] - grid -> z_vector[index_3];
                	if (delta_z > 0)
                		index_for_vertical_gradient = index_3 - NUMBER_OF_VECTORS_PER_LAYER;
                	else
                	{
                		if (layer_index == NUMBER_OF_LAYERS - 1)
                			index_for_vertical_gradient = index_3 - NUMBER_OF_VECTORS_PER_LAYER;
                		else
                			index_for_vertical_gradient = index_3 + NUMBER_OF_VECTORS_PER_LAYER;
                	}
                	vertical_gradient = (in_field[index_3] - in_field[index_for_vertical_gradient])/(grid -> z_vector[index_3] - grid -> z_vector[index_for_vertical_gradient]);
                	velocity_value_3 = in_field[index_3] + delta_z*vertical_gradient;
                	rhombus_circ = length_rescale_factor_0*grid -> normal_distance[index_0]*sign_0*velocity_value_0 + length_rescale_factor_1*grid -> normal_distance[index_1]*sign_1*velocity_value_1 + length_rescale_factor_2*grid -> normal_distance[index_2]*sign_2*velocity_value_2 + length_rescale_factor_3*grid -> normal_distance[index_3]*sign_3*velocity_value_3;
            	}
                else
                	rhombus_circ = grid -> normal_distance[index_0]*sign_0*in_field[index_0] + grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_2]*sign_2*in_field[index_2] + grid -> normal_distance[index_3]*sign_3*in_field[index_3];
                out_field[i] += 1.0/3*rhombus_circ/dualgrid -> area[edge_vector_index_dual_area];
            }
        }
        else
        {
            if (layer_index == 0 || layer_index == NUMBER_OF_LAYERS)
            {
                index_1 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_3 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 3];
                sign_1 = dualgrid -> h_curl_signs[4*h_index + 1];
                sign_3 = dualgrid -> h_curl_signs[4*h_index + 3];
                out_field[i] = 1/dualgrid -> area[i]*(grid -> normal_distance[index_1]*sign_1*in_field[index_1] + grid -> normal_distance[index_3]*sign_3*in_field[index_3]);
            }
            else
            {
                index_0 = layer_index*NUMBER_OF_VECTORS_PER_LAYER - NUMBER_OF_VECTORS_H + dualgrid -> h_curl_indices[4*h_index + 0];
                index_1 = layer_index*NUMBER_OF_VECTORS_PER_LAYER + dualgrid -> h_curl_indices[4*h_index + 1];
                index_2 = layer_index*NUMBER_OF_VECTORS_PER_LAYER - NUMBER_OF_VECTORS_H + dualgrid -> h_curl_indices[4*h_index + 2];
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
                if (layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
                {
                    retval = horizontal_covariant_normalized(in_field, layer_index, dualgrid -> h_curl_indices[4*h_index + 2], grid, &covar_0);
                    if (retval != 0)
                        printf("Error in horizontal_covariant_normalized called at position 0 from curl.\n");
                }
                covar_2 = in_field[index_2];
                if (layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS + 1)
                {
                    retval = horizontal_covariant_normalized(in_field, layer_index - 1, dualgrid -> h_curl_indices[4*h_index + 2], grid, &covar_2);
                    if (retval != 0)
                        printf("Error in horizontal_covariant_normalized called at position 1 from curl.\n");
                }
                out_field[i] = 1/dualgrid -> area[i]*(dist_0*sign_0*covar_0 + dist_1*sign_1*in_field[index_1] + dist_2*sign_2*covar_2 + dist_3*sign_3*in_field[index_3]);
            }
        }
    }
    return 0;
}












