/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This file contains functions that compute fundamental properties of the vertical grid.
*/

#include "include.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geos95.h"

int set_z_scalar(double z_scalar[], double z_surface[], int NO_OF_ORO_LAYERS, double TOA, double stretching_parameter, int VERT_GRID_TYPE)
{
	/*
	This function sets the z coordinates of the scalar data points.
	*/
	
	double z_vertical_vector_pre[NO_OF_LAYERS + 1];
	// the heights are defined according to z_k = A_k + B_k*z_surface with A_0 = TOA, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
	double A, B, sigma_z, z_rel, max_oro;
	// loop over all columns
    for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
    {
    	// filling up z_vertical_vector_pre
		for (int j = 0; j < NO_OF_LAYERS + 1; ++j)
		{
			z_rel = 1 - (j + 0.0)/NO_OF_LAYERS; // z/TOA
			sigma_z = pow(z_rel, stretching_parameter);
			A = sigma_z*TOA; // the height without orography
			// B corrects for orography
			if (j >= NO_OF_LAYERS - NO_OF_ORO_LAYERS && VERT_GRID_TYPE == 0)
			{
				B = (j - (NO_OF_LAYERS - NO_OF_ORO_LAYERS) + 0.0)/NO_OF_ORO_LAYERS;
			}
			else
			{
				B = 0;
			}
			z_vertical_vector_pre[j] = A + B*z_surface[h_index];
		}
		
		// doing a check
		if (h_index == 0 && VERT_GRID_TYPE == 0)
		{
			max_oro = z_surface[find_max_index(z_surface, NO_OF_SCALARS_H)];
			if (max_oro >= z_vertical_vector_pre[NO_OF_LAYERS - NO_OF_ORO_LAYERS])
			{
				printf("Maximum of orography larger or equal to the height of the lowest flat level.\n");
				printf("Aborting.\n");
				exit(1);
			}
		}
		
		// placing the scalar points in the middle between the preliminary values of the adjacent levels
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			z_scalar[layer_index*NO_OF_SCALARS_H + h_index] = 0.5*(z_vertical_vector_pre[layer_index] + z_vertical_vector_pre[layer_index + 1]);
    	}
    }
    return 0;
}

int set_scalar_shading_indices(double z_scalar[], double z_surface[], int no_of_shaded_points_scalar[])
{
	int counter;
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		counter = 0;
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (z_scalar[j*NO_OF_SCALARS_H + i] < z_surface[i])
			{
				++counter;
			}
		}
		no_of_shaded_points_scalar[i] = counter;
	}
	return 0;
}

int set_vector_shading_indices(int from_index[], int to_index[], int no_of_shaded_points_scalar[], int no_of_shaded_points_vector[])
{
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		no_of_shaded_points_vector[i] = fmax(no_of_shaded_points_scalar[from_index[i]], no_of_shaded_points_scalar[to_index[i]]);
	}
	return 0;
}

int set_z_vector_and_normal_distance(double z_vector[], double z_scalar[], double normal_distance[], double latitude_scalar[], double longitude_scalar[], int from_index[], int to_index[], double TOA, int VERT_GRID_TYPE, double z_surface[])
{
	/*
	calculates the vertical position of the vector points
	as well as the normal distances of the primal grid
	*/
	
	int layer_index, h_index, upper_index, lower_index;
	double *lowest_thicknesses = malloc(NO_OF_SCALARS_H*sizeof(double));
	#pragma omp parallel for private(layer_index, h_index, upper_index, lower_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        // horizontal grid points
        if (h_index >= NO_OF_SCALARS_H)
        {
        	// placing the vector vertically in the middle between the two adjacent scalar points
            z_vector[i]
            = 0.5*(z_scalar[layer_index*NO_OF_SCALARS_H + from_index[h_index - NO_OF_SCALARS_H]]
            + z_scalar[layer_index*NO_OF_SCALARS_H + to_index[h_index - NO_OF_SCALARS_H]]);
            // calculating the horizontal distance
            normal_distance[i]
            = calculate_distance_h(
            latitude_scalar[from_index[h_index - NO_OF_SCALARS_H]],
            longitude_scalar[from_index[h_index - NO_OF_SCALARS_H]],
            latitude_scalar[to_index[h_index - NO_OF_SCALARS_H]],
            longitude_scalar[to_index[h_index - NO_OF_SCALARS_H]],
            RADIUS + z_vector[i]);
        }
        else
        {
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            // highest level
            if (layer_index == 0)
			{
            	z_vector[i] = TOA;
                normal_distance[i] = TOA - z_scalar[lower_index];
			}
			// lowest level
            else if (layer_index == NO_OF_LAYERS)
			{
				if (VERT_GRID_TYPE == 0)
				{
					z_vector[i] = z_surface[h_index];
                }
				if (VERT_GRID_TYPE == 1)
				{
					z_vector[i] = 0;
                }
                normal_distance[i] = z_scalar[upper_index] - z_vector[i];
                lowest_thicknesses[h_index] = z_vector[i - NO_OF_VECTORS_PER_LAYER] - z_vector[i];
			}
			// inner levels
            else
			{
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
				// placing the vertical vector in the middle between the two adjacent scalar points
				z_vector[i] = z_scalar[lower_index] + 0.5*normal_distance[i];
			}
        }
    }
    
    // checks
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
		if (fabs(check_sum/(TOA - z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i]) - 1) > EPSILON_SECURITY)
		{
			printf("Problem 0 with vertical grid structure.\n");
			exit(1);
		}
	}
	double max_thick, min_thick, thick_rel;
	min_thick = lowest_thicknesses[find_min_index(lowest_thicknesses, NO_OF_SCALARS_H)];
	max_thick = z_vector[0] - z_vector[NO_OF_VECTORS_PER_LAYER];
	thick_rel = max_thick/min_thick;
	printf("ratio of maximum to minimum layer thickness (including orography): %lf\n", thick_rel);
	free(lowest_thicknesses);
    return 0;
}

int calculate_vertical_faces(double area[], double z_vector_dual[], double normal_distance_dual[], double TOA)
{
	/*
	This function calculates the vertical faces.
	*/
	
	int layer_index, h_index, dual_vector_index;
	double base_distance, radius_0, radius_1;
	#pragma omp parallel for private(layer_index, h_index, dual_vector_index, base_distance, radius_0, radius_1)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            dual_vector_index = (layer_index + 1)*NO_OF_DUAL_VECTORS_PER_LAYER + h_index - NO_OF_SCALARS_H;
            radius_0 = RADIUS + z_vector_dual[dual_vector_index];
            radius_1 = RADIUS + z_vector_dual[dual_vector_index - NO_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index];
            area[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    
    // checks
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        if (area[i] <= 0)
		{
            printf("It is area <= 0 at some point, position 0.\n");
			exit(1);
		}
    }
    double check_area, wished_result;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	wished_result = calculate_vertical_face(normal_distance_dual[NO_OF_LAYERS*NO_OF_DUAL_VECTORS_PER_LAYER + i], RADIUS + z_vector_dual[NO_OF_LAYERS*NO_OF_DUAL_VECTORS_PER_LAYER + i], RADIUS + TOA);
    	check_area = 0;
    	for (int j = 0; j < NO_OF_LAYERS; ++j)
    	{
    		check_area += area[NO_OF_SCALARS_H + i + j*NO_OF_VECTORS_PER_LAYER];
    	}
    	if(fabs(check_area/wished_result - 1) > EPSILON_SECURITY)
    	{
    		printf("Error with vertical faces. Coefficient which should be zero has value %lf.\n", check_area/wished_result - 1);
    		exit(1);
    	}
    }
    return 0;
}

int set_z_scalar_dual(double z_scalar_dual[], double z_vector[], int from_index[], int to_index[], int vorticity_indices_triangles[], double TOA)
{
	/*
	This function sets the z coordinates of the dual scalar points.
	*/
	
	int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_DUAL_SCALARS; ++i)
    {
		layer_index = i/NO_OF_DUAL_SCALARS_H;
		h_index = i - layer_index*NO_OF_DUAL_SCALARS_H;
		z_scalar_dual[i]
		= 1.0/6*(
		z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + from_index[vorticity_indices_triangles[3*h_index + 0]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + from_index[vorticity_indices_triangles[3*h_index + 1]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + from_index[vorticity_indices_triangles[3*h_index + 2]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + to_index[vorticity_indices_triangles[3*h_index + 0]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + to_index[vorticity_indices_triangles[3*h_index + 1]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + to_index[vorticity_indices_triangles[3*h_index + 2]]]);
    }
    
    // checks
	int min_oro_index = find_min_index(z_vector, NO_OF_VECTORS);
	double min_oro = z_vector[min_oro_index];
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

int set_volume(double volume[], double z_vector[], double area[], int from_index[], int to_index[], double TOA, int vorticity_indices_triangles[])
{
	/*
	This function computes the volumes of the grid boxes.
	*/
	
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
    
    // checks
    // 1.) grid box volumes always need to be positive
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        if (volume[i] <= 0)
		{
            printf("volume contains a non-positive value.\n");
			exit(1);
		}
    }
    
    // 2.) check if the sum of all grid box volumes is the same as the volume of the whole model atmosphere
    volume_sum_ideal = 0;
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
    	volume_sum_ideal += find_volume(area[NO_OF_VECTORS - NO_OF_SCALARS_H + i], RADIUS + z_vector[NO_OF_VECTORS- NO_OF_SCALARS_H + i], RADIUS + TOA);
    }
    if (fabs(volume_sum/volume_sum_ideal - 1) > EPSILON_SECURITY)
	{
        printf("Sum of volumes of grid boxes does not match volume of entire atmosphere.\n");
		exit(1);
	}
	return 0;
}

int set_area_dual(double area_dual[], double z_vector_dual[], double normal_distance[], double z_vector[], int from_index[], int to_index[], double triangle_face_unit_sphere[], double TOA)
{
	int layer_index, h_index, primal_vector_index;
	double radius_0, radius_1, base_distance;
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NO_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_H)
        {
            area_dual[i] = pow(RADIUS + z_vector_dual[i], 2)*triangle_face_unit_sphere[h_index - NO_OF_VECTORS_H];
        }
        else
        {
        	if (layer_index == 0)
        	{
		        primal_vector_index = NO_OF_SCALARS_H + h_index;
		        radius_0 = RADIUS + z_vector[primal_vector_index];
		        radius_1 = RADIUS + TOA;
		        base_distance = normal_distance[primal_vector_index];
        	}
        	else if (layer_index == NO_OF_LAYERS)
        	{
		        primal_vector_index = NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + h_index;
		        radius_0 = RADIUS + 0.5*(z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + from_index[h_index]] + z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + to_index[h_index]]);
		        radius_1 = RADIUS + z_vector[primal_vector_index];
		        base_distance = normal_distance[primal_vector_index]*radius_0/radius_1;
        	}
        	else
        	{
		        primal_vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
		        radius_0 = RADIUS + z_vector[primal_vector_index];
		        radius_1 = RADIUS + z_vector[primal_vector_index - NO_OF_VECTORS_PER_LAYER];
		        base_distance = normal_distance[primal_vector_index];
        	}
            area_dual[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        if (area_dual[i] <= 0)
		{
            printf("area_dual contains a non-positive value.\n");
			exit(1);
		}
    }
    double check_area, wished_result;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        radius_0 = RADIUS + 0.5*(z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + from_index[i]] + z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + to_index[i]]);
        primal_vector_index = NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + i;
        radius_1 = RADIUS + z_vector[primal_vector_index];
        base_distance = normal_distance[primal_vector_index]*radius_0/radius_1;
    	wished_result = calculate_vertical_face(base_distance, radius_0, RADIUS + TOA);
    	check_area = 0;
    	for (int j = 0; j < NO_OF_LEVELS; ++j)
    	{
    		check_area += area_dual[i + j*NO_OF_DUAL_VECTORS_PER_LAYER];
    	}
    	if(fabs(check_area/wished_result - 1) > EPSILON_SECURITY)
    	{
    		printf("Error with dual vertical faces. Coefficient which should be zero has value %lf.\n", check_area/wished_result - 1);
    		exit(1);
    	}
    }
    return 0;
}

int set_gravity_potential(double z_scalar[], double gravity_potential[], double GRAVITY_MEAN_SFC_ABS)
{
	/*
	Thi function computes the gravity potential.
	*/
	
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	gravity_potential[i] = -GRAVITY_MEAN_SFC_ABS*(RADIUS*RADIUS/(RADIUS + z_scalar[i]) - RADIUS);
    }
	return 0;
}

int map_hor_area_to_half_levels(double area[], double z_vector[], double pent_hex_face_unity_sphere[])
{
	/*
	This function maps the horizontal areas from the surface of the unit sphere to the half-levels.
	*/
	
	int layer_index, h_index;
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_SCALARS_H)
        {
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(RADIUS + z_vector[i], 2);
        }
    }
	return 0;
}

int calc_z_vector_dual_and_normal_distance_dual(double z_vector_dual[], double normal_distance_dual[], double z_scalar_dual[], double TOA, int from_index[],
int to_index[], double z_vector[], int from_index_dual[], int to_index_dual[], double latitude_scalar_dual[], double longitude_scalar_dual[], int vorticity_indices_triangles[])
{
	/*
	This function sets the z coordinates of the dual vector points as well as the normal distances of the dual grid.
	*/
	
	int layer_index, h_index, upper_index, lower_index;
	#pragma omp parallel for private(layer_index, h_index, upper_index, lower_index)
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NO_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_H)
        {
            upper_index = h_index - NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_SCALARS_H;
            lower_index = h_index - NO_OF_VECTORS_H + (layer_index + 1)*NO_OF_DUAL_SCALARS_H;
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
			z_vector_dual[i] = 1.0/3*(z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + vorticity_indices_triangles[3*(h_index - NO_OF_VECTORS_H) + 0]]
			+ z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + vorticity_indices_triangles[3*(h_index - NO_OF_VECTORS_H) + 1]]
			+ z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + vorticity_indices_triangles[3*(h_index - NO_OF_VECTORS_H) + 2]]);
        }
        else
        {
			if (layer_index == 0)
			{
				z_vector_dual[i] = TOA;
			}
			else if (layer_index == NO_OF_LAYERS)
			{
				z_vector_dual[i] = 0.5*(z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + from_index[h_index]] + z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + to_index[h_index]]);
			}
			else
			{
				z_vector_dual[i] = 0.5*(z_vector[NO_OF_SCALARS_H + h_index + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER] + z_vector[NO_OF_SCALARS_H + h_index + layer_index*NO_OF_VECTORS_PER_LAYER]);
			}
            normal_distance_dual[i] = calculate_distance_h(latitude_scalar_dual[from_index_dual[h_index]], longitude_scalar_dual[from_index_dual[h_index]],
            latitude_scalar_dual[to_index_dual[h_index]], longitude_scalar_dual[to_index_dual[h_index]],
            RADIUS + z_vector_dual[i]);
        }
    }
    
    // checks
	int index_vector_for_dual_scalar_z[3];
	double check_sum;
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        if (normal_distance_dual[i] <= 0)
		{
            printf("normal_distance_dual contains a non-positive value.\n");
            printf("Aborting.\n");
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
		find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_triangles, i, index_vector_for_dual_scalar_z);
		if (fabs(check_sum/(TOA
		- 1.0/3*(z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]]
		+ z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]]
		+ z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]])) - 1) > EPSILON_SECURITY)
		{
			printf("Problem 1 with vertical grid structure.\n");
			printf("Aborting.\n");
			exit(1);
		}
	}
	return 0;
}

int slopes(double z_scalar[], int from_index[], int to_index[], double normal_distance[], double slope[])
{
	/*
	calculates the slopes of the terrain following coordinates
	*/
	int layer_index, h_index;
	double delta_x, delta_z;
	#pragma omp parallel for private(layer_index, h_index, delta_x, delta_z)
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		// At the vertical vector points, the slopes are set to zero.
		if (h_index < NO_OF_SCALARS_H)
		{
			slope[i] = 0;
		}
		else
		{
			delta_x = normal_distance[layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
			delta_z
			= z_scalar[layer_index*NO_OF_SCALARS_H + to_index[h_index - NO_OF_SCALARS_H]]
			- z_scalar[layer_index*NO_OF_SCALARS_H + from_index[h_index - NO_OF_SCALARS_H]];
			slope[i] = delta_z/delta_x;
			// Checking for unrealistic values.
			if (fabs(slope[i]) > 1)
			{
				printf("Problem with slopes.\n");
				printf("Aborting.\n");
				exit(1);
			}
		}
	}
	return 0;
}














