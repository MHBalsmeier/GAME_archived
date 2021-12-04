/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This function is a collection of some helper functions that are needed for the grid generator.
*/

#include <netcdf.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <geos95.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "include.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int set_f_vec(double latitude_vector[], double direction[], double direction_dual[], double f_vec[])
{
	/*
	This function sets the Coriolis vector (vertical at horizontal primal vector points,
	horizontal at horizontal dual vector points).
	*/
	
	#pragma omp parallel for
    for (int i = 0; i < 2*NO_OF_VECTORS_H; ++i)
    {
    	// horizontal component at dual vector points
        if (i < NO_OF_VECTORS_H)
        {
        	f_vec[i] = 2*OMEGA*cos(latitude_vector[i])*sin(direction_dual[i]);
    	}
        // vertical component at primal vector points
        else if (i < 2*NO_OF_VECTORS_H)
        {
        	f_vec[i] = 2*OMEGA*sin(latitude_vector[i - NO_OF_VECTORS_H]);
    	}
    }
 	return 0;   
}

int set_orography(int RES_ID, int ORO_ID, double z_surface[])
{
	/*
	This function reads the orography from a netcdf file.
	*/

    if (ORO_ID != 0)
    {
    	int retval, z_surface_id, ncid;
		int ORO_FILE_LENGTH = 100;
		char *ORO_FILE_PRE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE_PRE, "../surface_generator/surface_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
		ORO_FILE_LENGTH = strlen(ORO_FILE_PRE);
		free(ORO_FILE_PRE);
		char *ORO_FILE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE, "../surface_generator/surface_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);	    
		if ((retval = nc_open(ORO_FILE, NC_NOWRITE, &ncid)))
		    ERR(retval);
		if ((retval = nc_inq_varid(ncid, "z_surface", &z_surface_id)))
		    ERR(retval);		
		if ((retval = nc_get_var_double(ncid, z_surface_id, &z_surface[0])))
		    ERR(retval);
		if ((retval = nc_close(ncid)))
		    ERR(retval);
    }
    else
    {
    	#pragma omp parallel for
    	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    	{
    		z_surface[i] = 0;
    	}
    }
	return 0;
}

int calc_vorticity_indices_triangles(int from_index_dual[], int to_index_dual[], double direction[], double direction_dual[], int vorticity_indices_triangles[], double ORTH_CRITERION_DEG, int vorticity_signs_pre[])
{
	int counter, sign;
	double direction_change;
	#pragma omp parallel for private(counter, sign, direction_change)
    for (int i = 0; i < NO_OF_DUAL_SCALARS_H; ++i)
    {
        counter = 0;
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            if (from_index_dual[j] == i || to_index_dual[j] == i)
            {
                vorticity_indices_triangles[3*i + counter] = j;
                sign = 1;
                if (from_index_dual[j] == i)
                {
                    direction_change = find_turn_angle(direction_dual[j], direction[j]);
                    if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                    {
                        sign = -1;
                    }
                }
                if (to_index_dual[j] == i)
                {
                    direction_change = find_turn_angle(direction_dual[j], direction[j]);
                    if (rad2deg(direction_change) > ORTH_CRITERION_DEG)
                    {
                        sign = -1;
                    }
                }
                vorticity_signs_pre[3*i + counter] = sign;
                ++counter;
            }
        }
        if (counter != 3)
		{
            printf("Trouble detected, place 0.\n");
			exit(1);
		}
    }
	return 0;
}

int write_statistics_file(double pent_hex_face_unity_sphere[], double normal_distance[], double normal_distance_dual[], char statistics_file_name[])
{
	/*
	This function writes out statistical properties of the grid to a text file.
	*/
    double area_max, area_min, normal_distance_h_min, normal_distance_h_max, normal_distance_dual_h_min, normal_distance_dual_h_max;
    area_min = pent_hex_face_unity_sphere[find_min_index(pent_hex_face_unity_sphere, NO_OF_SCALARS_H)];
    area_max = pent_hex_face_unity_sphere[find_max_index(pent_hex_face_unity_sphere, NO_OF_SCALARS_H)];
    double *horizontal_distance = malloc(NO_OF_VECTORS_H*sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	horizontal_distance[i] = normal_distance[NO_OF_SCALARS_H + i];
    }
    normal_distance_h_min = horizontal_distance[find_min_index(horizontal_distance, NO_OF_VECTORS_H)];
    normal_distance_h_max = horizontal_distance[find_max_index(horizontal_distance, NO_OF_VECTORS_H)];
    double *horizontal_distance_dual = malloc(NO_OF_VECTORS_H*sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	horizontal_distance_dual[i] = normal_distance_dual[i];
    }
    normal_distance_dual_h_min = horizontal_distance_dual[find_min_index(horizontal_distance_dual, NO_OF_VECTORS_H)];
    normal_distance_dual_h_max = horizontal_distance_dual[find_max_index(horizontal_distance_dual, NO_OF_VECTORS_H)];
    FILE *statistics_file = fopen(statistics_file_name, "w");
    fprintf(statistics_file, "Ratio of minimum to maximum area: %lf\n", area_min/area_max);
   	fprintf(statistics_file, "Shortest horizontal normal distance (highest layer): %lf m.\n", normal_distance_h_min);
    fprintf(statistics_file, "Longest horizontal normal distance (highest layer): %lf m.\n", normal_distance_h_max);
    fprintf(statistics_file, "Ratio of shortest to longest horizontal normal distance: %lf.\n", normal_distance_h_min/normal_distance_h_max);
    fprintf(statistics_file, "Shortest horizontal normal distance dual (highest level): %lf m.\n", normal_distance_dual_h_min);
    fprintf(statistics_file, "Longest horizontal normal distance dual (highest level): %lf m.\n", normal_distance_dual_h_max);
    fprintf(statistics_file, "Ratio of shortest to longest dual horizontal normal distance: %lf.\n", normal_distance_dual_h_min/normal_distance_dual_h_max);
    fclose(statistics_file);
    free(horizontal_distance);
    free(horizontal_distance_dual);
	return 0;
}












