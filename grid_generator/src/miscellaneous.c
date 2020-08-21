/*
This source file is part of the General Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <netcdf.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "grid_generator.h"
#include "enum.h"
#include "geos95.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int set_f_vec(double latitude_vector[], double direction[], double direction_dual[], double f_vec[])
{
	#pragma omp parallel for
    for (int i = 0; i < 3*NO_OF_VECTORS_H; ++i)
    {
    	// horizontal component at dual vector points
        if (i < NO_OF_VECTORS_H)
        	f_vec[i] = 2*OMEGA*cos(latitude_vector[i])*sin(direction_dual[i]);
        // vertical component at primal vector points
        else if (i < 2*NO_OF_VECTORS_H)
        	f_vec[i] = 2*OMEGA*sin(latitude_vector[i - NO_OF_VECTORS_H]);
    	// preparation of horizontal non-traditional Coriolis term
        else
        	f_vec[i] = 2*OMEGA*cos(latitude_vector[i - 2*NO_OF_VECTORS_H])*cos(direction[i - 2*NO_OF_VECTORS_H]);
    }
 	return 0;   
}

int find_angle_change(double angle_0, double angle_1, double *result)
{
    double result_pre = angle_1 - angle_0;
    if (result_pre > M_PI)
        result_pre = result_pre - 2*M_PI;
    if (result_pre < -M_PI)
        result_pre = result_pre + 2*M_PI;
    *result = result_pre;
    return 0;
}

int set_orography(int RES_ID, int ORO_ID, double z_surface[])
{
    if (ORO_ID != 0)
    {
    	int retval, z_surface_id, ncid;
		int ORO_FILE_LENGTH = 100;
		char *ORO_FILE_PRE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE_PRE, "../orography_generator/nc_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
		ORO_FILE_LENGTH = strlen(ORO_FILE_PRE);
		free(ORO_FILE_PRE);
		char *ORO_FILE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE, "../orography_generator/nc_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);	    
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
    	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    	{
    		z_surface[i] = 0;
    	}
    }
	return 0;
}

int check_for_orthogonality(double direction[], double direction_dual[], double ORTH_CRITERION_DEG)
{
	double direction_change;
	#pragma omp parallel for private(direction_change)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        find_angle_change(direction[i], direction_dual[i], &direction_change);
        if (fabs(rad2deg(direction_change)) < ORTH_CRITERION_DEG || fabs(rad2deg(direction_change)) > 90 + (90 - ORTH_CRITERION_DEG))
		{
            printf("Grid non-orthogonal: Intersection angle of %lf degrees detected.\n", fabs(rad2deg(direction_change)));
			exit(1);
		}
    }
    return 0;
}

int calc_vorticity_indices_pre(int from_index_dual[], int to_index_dual[], double direction[], double direction_dual[], int vorticity_indices_pre[], double ORTH_CRITERION_DEG, int vorticity_signs_pre[])
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
                vorticity_indices_pre[3*i + counter] = j;
                sign = 1;
                if (from_index_dual[j] == i)
                {
                    find_angle_change(direction_dual[j], direction[j], &direction_change);
                    if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                        sign = -1;
                }
                if (to_index_dual[j] == i)
                {
                    find_angle_change(direction_dual[j], direction[j], &direction_change);
                    if (rad2deg(direction_change) > ORTH_CRITERION_DEG)
                        sign = -1;
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
    double area_max, area_min, normal_distance_h_min, normal_distance_h_max, normal_distance_dual_h_min, normal_distance_dual_h_max;
    area_min = pent_hex_face_unity_sphere[find_min_index(pent_hex_face_unity_sphere, NO_OF_SCALARS_H)];
    area_max = pent_hex_face_unity_sphere[find_max_index(pent_hex_face_unity_sphere, NO_OF_SCALARS_H)];
    double *horizontal_distance = malloc(NO_OF_VECTORS_H*sizeof(double));
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    	horizontal_distance[i] = normal_distance[NO_OF_SCALARS_H + i];
    normal_distance_h_min = horizontal_distance[find_min_index(horizontal_distance, NO_OF_VECTORS_H)];
    normal_distance_h_max = horizontal_distance[find_max_index(horizontal_distance, NO_OF_VECTORS_H)];
    double *horizontal_distance_dual = malloc(NO_OF_VECTORS_H*sizeof(double));
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    	horizontal_distance_dual[i] = normal_distance_dual[i];
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












