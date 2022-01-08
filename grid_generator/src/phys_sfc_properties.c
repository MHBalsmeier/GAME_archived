/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
With this program, orographies can be produced.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include <geos95.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define P_0 100000.0

const double MOUNTAIN_HEIGHT = 10e3;
const double MOUNTAIN_FWHM = 1000e3;

int set_sfc_properties(double latitude_scalar[], double longitude_scalar[], double roughness_length[], double sfc_albedo[],
double sfc_rho_c[], double t_conductivity[], double oro[], int is_land[], int oro_id)
{
	double *oro_unfiltered = malloc(NO_OF_SCALARS_H*sizeof(double));
	
	int ncid, retval, is_land_id;
	if (oro_id == 2)
	{
		char is_land_file_pre[200];
		sprintf(is_land_file_pre, "real/B%d_is_land.nc", RES_ID);
		char is_land_file[strlen(is_land_file_pre) + 1];
		strcpy(is_land_file, is_land_file_pre);
		if ((retval = nc_open(is_land_file, NC_NOWRITE, &ncid)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "is_land", &is_land_id)))
			ERR(retval);
		if ((retval = nc_get_var_int(ncid, is_land_id, &is_land[0])))
			ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
	}
  
  	// reading the ETOPO orography
   	int lat_in_id, lon_in_id, z_in_id;
   	int no_of_lat_points = 10801;
   	int no_of_lon_points = 21601;
	double *latitude_input = malloc(no_of_lat_points*sizeof(double));
   	double *longitude_input = malloc(no_of_lon_points*sizeof(double));
	int (*z_input)[no_of_lon_points] = malloc(sizeof(int[no_of_lat_points][no_of_lon_points]));
	if (oro_id == 2)
	{
		if ((retval = nc_open("real/etopo.nc", NC_NOWRITE, &ncid)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "y", &lat_in_id)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "x", &lon_in_id)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "z", &z_in_id)))
			ERR(retval);
		if ((retval = nc_get_var_double(ncid, lat_in_id, &latitude_input[0])))
			ERR(retval);
		if ((retval = nc_get_var_double(ncid, lon_in_id, &longitude_input[0])))
			ERR(retval);
		if ((retval = nc_get_var_int(ncid, z_in_id, &z_input[0][0])))
			ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
	}
	
    // setting the unfiltered orography
    int lat_index, lon_index;
    double sigma_mountain = MOUNTAIN_FWHM/pow(8*log(2), 0.5); // only for oro_id == 1
    double distance;
	#pragma omp parallel for private(distance, lat_index, lon_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// default
		oro[i] = 0;
		oro_unfiltered[i] = 0;
		if (oro_id == 1)
		{
            distance = calculate_distance_h(latitude_scalar[i], longitude_scalar[i], 0, 0, RADIUS);
			oro[i] = MOUNTAIN_HEIGHT*exp(-pow(distance, 2)/(2*pow(sigma_mountain, 2)));
		}
		// real orography can only be different from zero at land points
		if (oro_id == 2 && is_land[i] == 1)
		{
	    	double *lat_distance_vector = malloc(no_of_lat_points*sizeof(double));
   			double *lon_distance_vector = malloc(no_of_lon_points*sizeof(double));
   			for (int j = 0; j < no_of_lat_points; ++j)
   			{
   				lat_distance_vector[j] = fabs(deg2rad(latitude_input[j]) - latitude_scalar[i]);
   			}
   			for (int j = 0; j < no_of_lon_points; ++j)
   			{
   				lon_distance_vector[j] = fabs(deg2rad(longitude_input[j]) - longitude_scalar[i]);
   			}
			lat_index = find_min_index(lat_distance_vector, no_of_lat_points);
			lon_index = find_min_index(lon_distance_vector, no_of_lon_points);
			oro_unfiltered[i] = z_input[lat_index][lon_index];
			
			// check
			if (oro_unfiltered[i] < -382 || oro_unfiltered[i] > 8850)
			{
				printf("Warning: value out of usual range.\n");		
			}
			
			// freeing the memory
			free(lat_distance_vector);
			free(lon_distance_vector);
		}
	}
	free(z_input);
	free(latitude_input);
	free(longitude_input);
	
	// smoothing the real orography
	int no_of_avg_points = 8;
	int min_indices_vector[no_of_avg_points];
	double distance_vector[NO_OF_SCALARS_H];
	if (oro_id == 2)
	{
		#pragma omp parallel for private(min_indices_vector, distance_vector)
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			// finding the distance to the other grid points
			for (int j = 0; j < NO_OF_SCALARS_H; ++j)
			{
				distance_vector[j] = calculate_distance_h(latitude_scalar[i], longitude_scalar[i], latitude_scalar[j], longitude_scalar[j], 1);
			}
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				min_indices_vector[j] = -1;
			}
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				min_indices_vector[j] = find_min_index_exclude(distance_vector, NO_OF_SCALARS_H, min_indices_vector, no_of_avg_points);
			}
			oro[i] = 0;
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				oro[i] += oro_unfiltered[min_indices_vector[j]]/no_of_avg_points;
			}
		}
	}
	free(oro_unfiltered);
	
	printf("minimum orography: %lf m\n", oro[find_min_index(oro, NO_OF_SCALARS_H)]);
	printf("maximum orography: %lf m\n", oro[find_max_index(oro, NO_OF_SCALARS_H)]);
	
	double c_p_water = 4184.0;
	double c_p_soil = 830.0;
	double albedo_water = 0.06;
	// setting the land surface albedo to 0.12 (compare Zdunkowski, Trautmann & Bott:
	// Radiation in the Atmosphere, 2007, p. 444)
	double albedo_soil = 0.12;
	double density_soil = 1442.0;
	double t_conductivity_soil = 7.5e-7;
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// ocean
		sfc_albedo[i] = albedo_water;
		sfc_rho_c[i] = RHO_WATER*c_p_water;
		roughness_length[i] = 0.08;
		// for water this is set to some land-typical value, will not be used anyway
		t_conductivity[i] = t_conductivity_soil;
		// land
		if (is_land[i] == 1)
		{
			sfc_albedo[i] = albedo_soil;
			sfc_rho_c[i] = density_soil*c_p_soil;
			roughness_length[i] = 0.2;
		}
	}
	
	return 0;
}











