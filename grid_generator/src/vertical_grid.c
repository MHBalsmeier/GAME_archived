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
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		for (int j = 0; j < NO_OF_LAYERS + 1; ++j)
		{
			if (j >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				z_vertical_vector_pre[j] = z_surface[h_index] + (z_oro_off - z_surface[h_index])/NO_OF_ORO_LAYERS*(NO_OF_LAYERS - j);
			else
				z_vertical_vector_pre[j] = TOA - (TOA - z_oro_off)/(NO_OF_LAYERS - NO_OF_ORO_LAYERS)*j;
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
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_V)
        {
            dual_vector_index = NO_OF_LAYERS*NO_OF_DUAL_VECTORS_PER_LAYER + h_index - NO_OF_VECTORS_V;
            radius_1 = RADIUS + z_vector_dual[h_index - NO_OF_VECTORS_V + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER];
            radius_0 = RADIUS + z_vector_dual[h_index - NO_OF_VECTORS_V + (layer_index + 1)*NO_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index]/(RADIUS + z_vector_dual[dual_vector_index])*radius_0;
            area[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
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
	int layer_index, h_index, face_index, on_face_index, coord_0, coord_1, triangle_on_face_index, dual_scalar_index, dual_scalar_h_index, coord_0_points_amount;
	int index_vector_for_dual_scalar_z[3];
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_V)
        {
            if (h_index - NO_OF_VECTORS_V >= NO_OF_EDGES*(POINTS_PER_EDGE + 1))
            {
                face_index = (h_index - NO_OF_VECTORS_V - NO_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
                on_face_index = h_index - NO_OF_VECTORS_V - (NO_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
                triangle_on_face_index = on_face_index/3;
                find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
				dual_scalar_h_index = face_index*TRIANGLES_PER_FACE + 1 + 2*triangle_on_face_index + coord_1;
                dual_scalar_index = layer_index*NO_OF_DUAL_SCALARS_H + dual_scalar_h_index;
				find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index, index_vector_for_dual_scalar_z);
                z_scalar_dual[dual_scalar_index] = 1.0/3*(z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
				find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index - 1, index_vector_for_dual_scalar_z);
                z_scalar_dual[dual_scalar_index - 1] = 1.0/3*(z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                if (layer_index == NO_OF_LAYERS - 1)
                {
					find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + NO_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
					find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index - 1, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + NO_OF_DUAL_SCALARS_H - 1] = 1.0/3*(z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                }
                if (coord_0 == coord_0_points_amount - 1)
                {
					find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index + 1, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + 1] = 1.0/3*(z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    if (layer_index == NO_OF_LAYERS - 1)
            			z_scalar_dual[dual_scalar_index + 1 + NO_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    if (coord_1 == POINTS_PER_EDGE - 1)
                    {
						find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index + 2, index_vector_for_dual_scalar_z);
                		z_scalar_dual[dual_scalar_index + 2] = 1.0/3*(z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                        if (layer_index == NO_OF_LAYERS - 1)
            				z_scalar_dual[dual_scalar_index + 2 + NO_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    }
                }
            }
        }
    }
	int min_oro_index = find_min_index(z_surface, NO_OF_SCALARS_H);
	double min_oro = z_surface[min_oro_index];
	for (int i = 0; i < NO_OF_DUAL_SCALARS; ++i)
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
	return 0;
}

int set_volume(double volume[], double z_scalar_dual[], double z_vector[], double area[], int from_index[], int to_index[], double TOA, double z_surface[], int vorticity_indices_pre[])
{
	set_z_scalar_dual(z_scalar_dual, z_vector, from_index, to_index, vorticity_indices_pre, z_surface, TOA);
    double volume_sum, volume_sum_ideal, radius_0, radius_1, base_area;
    volume_sum = 0;
    int layer_index, h_index;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        radius_0 = RADIUS + z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        radius_1 = RADIUS + z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_0, radius_1);
        volume_sum += volume[i];
    }
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        if (volume[i] <= 0)
		{
            printf("volume contains a non-positive value.\n");
			exit(1);
		}
    }
    volume_sum_ideal = 0;
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    	volume_sum_ideal += find_volume(area[NO_OF_VECTORS - NO_OF_VECTORS_V + i], RADIUS + z_vector[NO_OF_VECTORS- NO_OF_VECTORS_V + i], RADIUS + TOA);
    if (fabs(volume_sum/volume_sum_ideal - 1) > 1e-12)
	{
        printf("Sum of volumes of grid boxes does not match volume of entire atmosphere.\n");
		exit(1);
	}
	return 0;
}

int calc_vertical_contravar_unit(double latitude_scalar[], double longitude_scalar[], int adjacent_vector_indices_h[], double vertical_contravar_unit[], double direction_dual[], double normal_distance_dual[], int from_index_dual[], int to_index_dual[], double z_scalar_dual[])
{
	int NO_OF_edges, layer_index, dual_vector_index;
    double x_value, y_value, z_value, normal_vector_global[3], local_basis_vector[3], local_x, local_y, local_z, tilting_angle, delta_x, delta_z, abs_value;
    for (int i = 0; i < NO_OF_ORO_LAYERS + 1; ++i)
    {
    	for (int j = 0; j < NO_OF_VECTORS_V; ++j)
    	{
    		if (j < NO_OF_PENTAGONS)
    		{
    			NO_OF_edges = 5;
			}
    		else
    		{
    			NO_OF_edges = 6;
			}
			normal_vector_global[0] = 0;
			normal_vector_global[1] = 0;
			normal_vector_global[2] = 0;
    		for (int k = 0; k < NO_OF_edges; ++k)
    		{
    			layer_index = i + NO_OF_LAYERS - NO_OF_ORO_LAYERS;
    			dual_vector_index = layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*j + k];
    			delta_x = normal_distance_dual[dual_vector_index];
    			delta_z = z_scalar_dual[layer_index*NO_OF_DUAL_SCALARS_H + to_index_dual[adjacent_vector_indices_h[6*j + k]]] - z_scalar_dual[layer_index*NO_OF_DUAL_SCALARS_H + from_index_dual[adjacent_vector_indices_h[6*j + k]]];
    			tilting_angle = atan(delta_z/delta_x);
    			local_x = -sin(tilting_angle)*cos(direction_dual[adjacent_vector_indices_h[6*j + k]]);
    			local_y = -sin(tilting_angle)*sin(direction_dual[adjacent_vector_indices_h[6*j + k]]);
    			local_z = cos(tilting_angle);	
    			abs_value = sqrt(pow(local_x, 2) + pow(local_y, 2) + pow(local_z, 2));
    			local_x = local_x/abs_value;
    			local_y = local_y/abs_value;
    			local_z = local_z/abs_value;
    			calc_local_i(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    			normal_vector_global[0] += local_x*local_basis_vector[0];
    			normal_vector_global[1] += local_x*local_basis_vector[1];
    			normal_vector_global[2] += local_x*local_basis_vector[2];
    			calc_local_j(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    			normal_vector_global[0] += local_y*local_basis_vector[0];
    			normal_vector_global[1] += local_y*local_basis_vector[1];
    			normal_vector_global[2] += local_y*local_basis_vector[2];
    			calc_local_k(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    			normal_vector_global[0] += local_z*local_basis_vector[0];
    			normal_vector_global[1] += local_z*local_basis_vector[1];
    			normal_vector_global[2] += local_z*local_basis_vector[2];
    		}
    		abs_value = sqrt(scalar_product_elementary(normal_vector_global, normal_vector_global));
    		normal_vector_global[0] = normal_vector_global[0]/abs_value;
    		normal_vector_global[1] = normal_vector_global[1]/abs_value;
    		normal_vector_global[2] = normal_vector_global[2]/abs_value;
    		calc_local_i(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    		x_value = scalar_product_elementary(normal_vector_global, local_basis_vector);
    		calc_local_j(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    		y_value = scalar_product_elementary(normal_vector_global, local_basis_vector);
    		calc_local_k(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    		z_value = scalar_product_elementary(normal_vector_global, local_basis_vector);
    		for (int k = 0; k < 3; ++k)
    		{
    			if (k == 0)
    			{
    				vertical_contravar_unit[i*3*NO_OF_VECTORS_V + 3*j + k] = x_value;
				}
    			if (k == 1)
    			{
    				vertical_contravar_unit[i*3*NO_OF_VECTORS_V + 3*j + k] = y_value;
				}
    			if (k == 2)
    			{
    				vertical_contravar_unit[i*3*NO_OF_VECTORS_V + 3*j + k] = z_value;
				}
				if (fabs(vertical_contravar_unit[i*3*NO_OF_VECTORS_V + 3*j + k]) > 1.000001)
				{
				    printf("fabs(vertical_contravar_unit) > 1 at some point.\n");
					exit(1);
				}
			}
    	}
    }
    return 0;
}

int rescale_area(double area[], double vertical_contravar_unit[])
{
    double area_rescale_factor, dz, dh, tilting_angle;
	double area_rescale_factor_warning_begin = 0.1;
	int contravar_unit_vector_index, layer_index, h_index;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	if (h_index < NO_OF_VECTORS_V && layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
    	{
    		contravar_unit_vector_index = h_index + (layer_index - (NO_OF_LAYERS - NO_OF_ORO_LAYERS))*NO_OF_VECTORS_V;
    		dz = vertical_contravar_unit[3*contravar_unit_vector_index + 2];
    		if (dz <= 0)
    			printf("Error in horizontal face orography rescale calculation, position 0.\n");
    		dh = sqrt(pow(vertical_contravar_unit[3*contravar_unit_vector_index + 0], 2) + pow(vertical_contravar_unit[3*contravar_unit_vector_index + 1], 2));
    		tilting_angle = atan(dh/dz);
    		area_rescale_factor = sqrt(1 + pow(tan(tilting_angle), 2));
    		area[i] = area_rescale_factor*area[i];
    		if (fabs(area_rescale_factor - 1) > area_rescale_factor_warning_begin)
    			printf("It is area_rescale_factor > area_rescale_factor_warning_begin at some point.\n");
    		if (area[i] <= 0)
    		{
    			printf("It is area <= 0 at some point, position 1.\n");
    			exit(1);
    		}
    	}
    }
	return 0;
}

int rescale_area_dual(double area_dual[], double z_vector[], int from_index_dual[], int to_index_dual[], double area_dual_pre[], double z_vector_dual[])
{
	int layer_index, h_index, primal_vector_index, dual_vector_index;
    double area_0, area_1, area_rescale_factor, area_ratio;
    for (int i = 0; i < NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS; ++i)
    {
    	layer_index = i/(2*NO_OF_VECTORS_H);
    	h_index = i - layer_index*2*NO_OF_VECTORS_H;
    	if (h_index < NO_OF_VECTORS_H)
    		area_dual[i] = area_dual_pre[layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + h_index];
    	else
    	{
    		primal_vector_index = NO_OF_VECTORS_V + h_index - NO_OF_VECTORS_H + layer_index*NO_OF_VECTORS_PER_LAYER;
    		dual_vector_index = NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[h_index - NO_OF_VECTORS_H];
    		area_0 = area_dual_pre[dual_vector_index];
    		area_rescale_factor = pow((RADIUS + z_vector[primal_vector_index])/(RADIUS + z_vector_dual[dual_vector_index]), 2);
    		area_0 = area_rescale_factor*area_0;
    		dual_vector_index = NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[h_index - NO_OF_VECTORS_H];
    		area_1 = area_dual_pre[dual_vector_index];
    		area_rescale_factor = pow((RADIUS + z_vector[primal_vector_index])/(RADIUS + z_vector_dual[dual_vector_index]), 2);
    		area_1 = area_rescale_factor*area_1;
    		area_ratio = area_0/area_1;
    		if (fabs(area_ratio - 1) > 0.3)
    		{
    			printf("Unrealistic area_ratio in rhombus area calculation, position 0.\n");
    			exit(1);
			}
    		area_dual[i] = area_0 + area_1;
    	}
		if (area_dual[i] <= 0)
		{
			printf("It is area_dual <= 0 at some point.\n");
			exit(1);
		}
    }
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	area_ratio = area_dual[i]/area_dual[i + 2*NO_OF_VECTORS_H];
    	if (fabs(area_ratio - 0.5) > 0.0001)
    	{
    		printf("Unrealistic value in area_ratio of area_dual, position 0.\n");
    		exit(1);
    	}
    }
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	area_ratio = area_dual[i + (NO_OF_LAYERS - 1)*2*NO_OF_VECTORS_H]/area_dual[i + (NO_OF_LAYERS - 1)*2*NO_OF_VECTORS_H + 2*NO_OF_VECTORS_H];
    	if (fabs(area_ratio - 2) > 0.001)
    	{
    		printf("Unrealistic value in area_ratio of area_dual, position 1.\n");
    		exit(1);
    	}
    }
    double check_area, area_sum, mean_z;
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
    	mean_z = 0;
    	for (int j = 0; j < NO_OF_VECTORS_H; ++j)
    		mean_z += z_vector[NO_OF_VECTORS_V + i*NO_OF_VECTORS_PER_LAYER + j];
    	mean_z = mean_z/NO_OF_VECTORS_H;
    	area_sum = 0;
    	for (int j = 0; j < NO_OF_VECTORS_H; ++j)
    		area_sum += 1.0/3*area_dual[NO_OF_VECTORS_H + i*2*NO_OF_VECTORS_H + j];
    	check_area = 4*M_PI*pow(RADIUS + mean_z, 2);
    	area_ratio = check_area/area_sum;
    	if (fabs(area_ratio - 1) > 0.00001)
    	{
    		printf("Unrealistic area_ratio in rhombus area calculation, position 1.\n");
    		exit(1);
    	}
    }
	return 0;
}

int set_gravity_potential(double z_scalar[], double gravity_potential[], double GRAVITY_MEAN_SFC_ABS)
{	
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	gravity_potential[i] = -GRAVITY_MEAN_SFC_ABS*(RADIUS*RADIUS/(RADIUS + z_scalar[i]) - RADIUS);
    }
	return 0;
}

int set_z_vector_and_normal_distance(double z_vector[], double z_surface[], double z_scalar[], double normal_distance[], double latitude_scalar[], double longitude_scalar[], int from_index[], int to_index[], double TOA)
{
	int layer_index, h_index, upper_index, lower_index;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_V)
        {
            z_vector[i] = 0.5*(z_scalar[layer_index*NO_OF_SCALARS_H + to_index[h_index - NO_OF_VECTORS_V]] + z_scalar[layer_index*NO_OF_SCALARS_H + from_index[h_index - NO_OF_VECTORS_V]]);
            if (z_vector[i] <= 0)
			{
                printf("z_vector contains a non-positive value at a horizontal grid point.\n");
				exit(1);
			}
            normal_distance[i] = calculate_distance_h(latitude_scalar[from_index[h_index - NO_OF_VECTORS_V]], longitude_scalar[from_index[h_index - NO_OF_VECTORS_V]], latitude_scalar[to_index[h_index - NO_OF_VECTORS_V]], longitude_scalar[to_index[h_index - NO_OF_VECTORS_V]], RADIUS + z_vector[i]);
        }
        else
        {
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            if (layer_index == 0)
			{
            	z_vector[i] = TOA;
                normal_distance[i] = TOA - z_scalar[lower_index];
			}
            else if (layer_index == NO_OF_LAYERS)
			{
				z_vector[i] = z_surface[h_index];
                normal_distance[i] = z_scalar[upper_index] - z_surface[h_index];
			}
            else
			{
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
				z_vector[i] = z_scalar[lower_index] + 0.5*normal_distance[i];
			}
            if (z_vector[i] < z_surface[h_index])
			{
                printf("z_vector lays below surface.\n");
				exit(1);
			}
        }
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        if (normal_distance[i] <= 0)
		{
            printf("normal_distance contains a non-positive value.\n");
			exit(1);
		}
    }
	double check_sum;
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		check_sum = 0;
		for (int j = 0; j < NO_OF_LEVELS; ++j)
		{
			check_sum += normal_distance[i + j*NO_OF_VECTORS_PER_LAYER];
		}
		if (fabs(check_sum/(TOA - z_surface[i]) - 1) > 1e-15)
		{
			printf("Problem 0 with vertical grid structure.\n");
			exit(1);
		}
	}
    return 0;
}

int map_area_to_sphere(double area[], double z_vector[], double pent_hex_face_unity_sphere[])
{
	int layer_index, h_index;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_VECTORS_V)
        {
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(RADIUS + z_vector[i], 2);
        }
    }
	return 0;
}

int calc_z_vector_dual_and_normal_distance_dual(double z_vector_dual[], double normal_distance_dual[], double z_scalar_dual[], double TOA, double z_surface[], int from_index[], int to_index[], double z_vector[], int from_index_dual[], int to_index_dual[], double latitude_scalar_dual[], double longitude_scalar_dual[], int vorticity_indices_pre[])
{
	int layer_index, h_index, upper_index, lower_index;
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NO_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_H)
        {
            upper_index = h_index - NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_SCALARS_H;
            lower_index = h_index - NO_OF_VECTORS_H + (layer_index + 1)*NO_OF_DUAL_SCALARS_H;
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
			z_vector_dual[i] = z_scalar_dual[lower_index] + 0.5*normal_distance_dual[i];
        }
        else
        {
			if (layer_index == 0)
				z_vector_dual[i] = TOA;
			else if (layer_index == NO_OF_LAYERS)
				z_vector_dual[i] = 0.5*(z_surface[from_index[h_index]] + z_surface[to_index[h_index]]);
			else
				z_vector_dual[i] = 0.5*(z_vector[NO_OF_VECTORS_V + h_index + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER] + z_vector[NO_OF_VECTORS_V + h_index + layer_index*NO_OF_VECTORS_PER_LAYER]);
            normal_distance_dual[i] = calculate_distance_h(latitude_scalar_dual[from_index_dual[h_index]], longitude_scalar_dual[from_index_dual[h_index]], latitude_scalar_dual[to_index_dual[h_index]], longitude_scalar_dual[to_index_dual[h_index]], RADIUS + z_vector_dual[i]);
        }
    }
	int index_vector_for_dual_scalar_z[3];
	double check_sum;
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        if (normal_distance_dual[i] <= 0)
		{
            printf("normal_distance_dual contains a non-positive value.\n");
			exit(1);
		}
	}
	for (int i = 0; i < NO_OF_DUAL_SCALARS_H; ++i)
	{
		check_sum = 0;
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			check_sum += normal_distance_dual[i + NO_OF_VECTORS_H + j*NO_OF_DUAL_VECTORS_PER_LAYER];
		}
		find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, i, index_vector_for_dual_scalar_z);
		if (fabs(check_sum/(TOA - 1.0/3*(z_surface[index_vector_for_dual_scalar_z[0]] + z_surface[index_vector_for_dual_scalar_z[1]] + z_surface[index_vector_for_dual_scalar_z[2]])) - 1) > 1e-15)
		{
			printf("Problem 1 with vertical grid structure.\n");
			exit(1);
		}
	}
	return 0;
}

int calc_area_dual_pre(double area_dual_pre[], double z_vector_dual[], double normal_distance[], double z_vector[], double triangle_face_unit_sphere[])
{
	int layer_index, h_index, primal_vector_index;
	double radius_0, radius_1, base_distance;
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NO_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_H)
        {
            area_dual_pre[i] = pow(RADIUS + z_vector_dual[i], 2)*triangle_face_unit_sphere[h_index - NO_OF_VECTORS_H];
        }
        else
        {
            if (layer_index == 0)
                radius_1 = RADIUS + z_vector_dual[i];
            else
                radius_1 = RADIUS + z_vector[NO_OF_VECTORS_V + h_index + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER];
            if (layer_index == NO_OF_LAYERS)
                radius_0 = RADIUS + z_vector_dual[i];
            else
                radius_0 = RADIUS + z_vector[NO_OF_VECTORS_V + h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
            primal_vector_index = (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + NO_OF_VECTORS_V + h_index;
            base_distance = normal_distance[primal_vector_index]/(RADIUS + z_vector[primal_vector_index])*radius_0;
            area_dual_pre[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    free(triangle_face_unit_sphere);
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        if (area_dual_pre[i] <= 0)
		{
            printf("area_dual_pre contains a non-positive value.\n");
			exit(1);
		}
    }
    return 0;
}






