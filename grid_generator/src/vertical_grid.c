/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int determine_z_scalar(double z_scalar[], double z_vertical_vector_pre[], double z_surface[], double z_oro_off, double TOA)
{
	int layer_index, h_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
		h_index = i - layer_index*NUMBER_OF_SCALARS_H;
		for (int j = 0; j < NUMBER_OF_LAYERS + 1; ++j)
		{
			if (j >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
				z_vertical_vector_pre[j] = z_surface[h_index] + (z_oro_off - z_surface[h_index])/NUMBER_OF_ORO_LAYERS*(NUMBER_OF_LAYERS - j);
			else
				z_vertical_vector_pre[j] = TOA - (TOA - z_oro_off)/(NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)*j;
		}
		z_scalar[i] = 0.5*(z_vertical_vector_pre[layer_index] + z_vertical_vector_pre[layer_index + 1]);
        if (z_scalar[i] <= 0)
		{
            printf("z_scalar contains a non-positive value.\n");
			exit(1);
		}
    }
    return 0;
}

int calculate_vertical_faces(double area[], double z_vector_dual[], double normal_distance_dual[])
{
	int layer_index, h_index, dual_vector_index;
	double base_distance, radius_0, radius_1;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            dual_vector_index = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_PER_LAYER + h_index - NUMBER_OF_VECTORS_V;
            radius_1 = RADIUS + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            radius_0 = RADIUS + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + (layer_index + 1)*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index]/(RADIUS + z_vector_dual[dual_vector_index])*radius_0;
            area[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        if (area[i] <= 0)
		{
            printf("It is area <= 0 at some point, position 0.\n");
			exit(1);
		}
    }
    return 0;
}

int set_z_scalar_dual(double z_scalar_dual[], double z_vector[], int from_index[], int to_index[], int vorticity_indices_pre[], double z_surface[], double TOA)
{
	int layer_index, h_index, face_index, on_face_index, retval, coord_0, coord_1, triangle_on_face_index, dual_scalar_index, dual_scalar_h_index, coord_0_points_amount;
	retval = 0;
	int index_vector_for_dual_scalar_z[3];
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            if (h_index - NUMBER_OF_VECTORS_V >= NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
            {
                face_index = (h_index - NUMBER_OF_VECTORS_V - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
                on_face_index = h_index - NUMBER_OF_VECTORS_V - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
                triangle_on_face_index = on_face_index/3;
                retval = find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
				dual_scalar_h_index = face_index*TRIANGLES_PER_FACE + 1 + 2*triangle_on_face_index + coord_1;
                dual_scalar_index = layer_index*NUMBER_OF_DUAL_SCALARS_H + dual_scalar_h_index;
				retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index, index_vector_for_dual_scalar_z);
                z_scalar_dual[dual_scalar_index] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
				retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index - 1, index_vector_for_dual_scalar_z);
                z_scalar_dual[dual_scalar_index - 1] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                if (layer_index == NUMBER_OF_LAYERS - 1)
                {
					retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + NUMBER_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
					retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index - 1, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + NUMBER_OF_DUAL_SCALARS_H - 1] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                }
                if (coord_0 == coord_0_points_amount - 1)
                {
					retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index + 1, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + 1] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    if (layer_index == NUMBER_OF_LAYERS - 1)
            			z_scalar_dual[dual_scalar_index + 1 + NUMBER_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    if (coord_1 == POINTS_PER_EDGE - 1)
                    {
						retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index + 2, index_vector_for_dual_scalar_z);
                		z_scalar_dual[dual_scalar_index + 2] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                        if (layer_index == NUMBER_OF_LAYERS - 1)
            				z_scalar_dual[dual_scalar_index + 2 + NUMBER_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    }
                }
            }
        }
    }
	int min_oro_index = find_min_index(z_surface, NUMBER_OF_SCALARS_H);
	double min_oro = z_surface[min_oro_index];
	for (int i = 0; i < NUMBER_OF_DUAL_SCALARS; ++i)
	{
		if (z_scalar_dual[i] > TOA)
		{
			printf("z_scalar_dual has a point above top of atmosphere.\n");
			exit(1);
		}
		if (z_scalar_dual[i] < min_oro)
		{
			printf("z_scalar_dual has a point below minimum of orography.\n");
			exit(1);
		}
	}
	return retval;
}

int set_volume(double volume[], double z_scalar_dual[], double z_vector[], double area[], int from_index[], int to_index[], double TOA, double z_surface[], int vorticity_indices_pre[])
{
	int retval = set_z_scalar_dual(z_scalar_dual, z_vector, from_index, to_index, vorticity_indices_pre, z_surface, TOA);
    double volume_sum, volume_sum_ideal, radius_0, radius_1, base_area;
    volume_sum = 0;
    int layer_index, h_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_0 = RADIUS + z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_1 = RADIUS + z_vector[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_0, radius_1);
        volume_sum += volume[i];
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (volume[i] <= 0)
		{
            printf("volume contains a non-positive value.\n");
			exit(1);
		}
    }
    volume_sum_ideal = 0;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    	volume_sum_ideal += find_volume(area[NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V + i], RADIUS + z_vector[NUMBER_OF_VECTORS- NUMBER_OF_VECTORS_V + i], RADIUS + TOA);
    if (fabs(volume_sum/volume_sum_ideal - 1) > 1e-12)
	{
        printf("Sum of volumes of grid boxes does not match volume of entire atmosphere.\n");
		exit(1);
	}
	return retval;
}



















