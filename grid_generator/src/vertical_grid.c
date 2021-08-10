/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This file contains functions that compute properties of the vertical grid.
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
	#pragma omp parallel for private(counter)
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
	#pragma omp parallel for
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
	double max_thick, min_thick, thick_rel;
	min_thick = lowest_thicknesses[find_min_index(lowest_thicknesses, NO_OF_SCALARS_H)];
	max_thick = z_vector[0] - z_vector[NO_OF_VECTORS_PER_LAYER];
	thick_rel = max_thick/min_thick;
	printf("ratio of maximum to minimum layer thickness (including orography): %lf\n", thick_rel);
	free(lowest_thicknesses);
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
	return 0;
}

int set_volume(double volume[], double z_vector[], double area[], int from_index[], int to_index[], double TOA, int vorticity_indices_triangles[])
{
	/*
	This function computes the volumes of the grid boxes.
	*/
	
    int layer_index, h_index;
    double radius_0, radius_1, base_area;
    #pragma omp parallel for private(layer_index, h_index, radius_0, radius_1, base_area)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        radius_0 = RADIUS + z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        radius_1 = RADIUS + z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_0, radius_1);
    }
	return 0;
}

int set_area_dual(double area_dual[], double z_vector_dual[], double normal_distance[], double z_vector[], int from_index[], int to_index[], double triangle_face_unit_sphere[], double TOA)
{
	int layer_index, h_index, primal_vector_index;
	double radius_0, radius_1, base_distance;
	#pragma omp parallel for private(layer_index, h_index, primal_vector_index, radius_0, radius_1, base_distance)
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
    return 0;
}

int set_area(double area[], double z_vector[], double z_vector_dual[], double normal_distance_dual[], double pent_hex_face_unity_sphere[])
{
	/*
	This function sets the areas of the grid boxes.
	*/
	
	int layer_index, h_index, dual_vector_index;
	double base_distance, radius_0, radius_1;
	#pragma omp parallel for private(layer_index, h_index, dual_vector_index, base_distance, radius_0, radius_1)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_SCALARS_H)
        {
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(RADIUS + z_vector[i], 2);
        }
        else
        {
            dual_vector_index = (layer_index + 1)*NO_OF_DUAL_VECTORS_PER_LAYER + h_index - NO_OF_SCALARS_H;
            radius_0 = RADIUS + z_vector_dual[dual_vector_index];
            radius_1 = RADIUS + z_vector_dual[dual_vector_index - NO_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index];
            area[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
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
	return 0;
}

int set_background_state(double z_scalar[], double gravity_potential[], double exner_bg[], double theta_bg[])
{
	/*
	This sets the hydrostatic background state.
	*/
	
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		
	}
	
	return 0;
}












