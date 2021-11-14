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
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define P_0 100000.0

const double MOUNTAIN_HEIGHT = 10e3;
const double MOUNTAIN_FWHM = 1000e3;

int main(int argc, char *argv[])
{
	int ORO_ID;
   	ORO_ID = strtod(argv[1], NULL);
   	if (ORO_ID < 1 || ORO_ID > 2)
   	{
   		printf("Error: oro_id must not be smaller than 1 or larger than 2.\n");
   		exit(1);
	}
   	char OUTPUT_FILE_PRE[200];
	sprintf(OUTPUT_FILE_PRE, "surface_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
   	char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
	strcpy(OUTPUT_FILE, OUTPUT_FILE_PRE);
	int ncid, scalar_h_dimid, oro_id, latitude_scalar_id, longitude_scalar_id;
	double *oro_unfiltered = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *oro = malloc(NO_OF_SCALARS_H*sizeof(double));
	char GEO_PROP_FILE_PRE[200];
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/grids/B%dL26T41152_O0_OL23_SCVT.nc", RES_ID);
	char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	int retval;
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double distance;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &latitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &longitude_scalar[0])))
        ERR(retval);
	if ((retval = nc_close(ncid)))
	  ERR(retval);
	
	// reading the land mask
	int *is_land = calloc(NO_OF_SCALARS_H, sizeof(int));
	int is_land_id;
	if (ORO_ID == 2)
	{
		char IS_LAND_FILE_PRE[200];
		sprintf(IS_LAND_FILE_PRE, "real/B%d_is_land.nc", RES_ID);
		char IS_LAND_FILE[strlen(IS_LAND_FILE_PRE) + 1];
		strcpy(IS_LAND_FILE, IS_LAND_FILE_PRE);
		if ((retval = nc_open(IS_LAND_FILE, NC_NOWRITE, &ncid)))
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
	if (ORO_ID == 2)
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
    double sigma_mountain = MOUNTAIN_FWHM/pow(8*log(2), 0.5); // only for ORO_ID == 1
	#pragma omp parallel for private(distance, lat_index, lon_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// default
		oro[i] = 0;
		oro_unfiltered[i] = 0;
		if (ORO_ID == 1)
		{
            distance = calculate_distance_h(latitude_scalar[i], longitude_scalar[i], 0, 0, RADIUS);
			oro[i] = MOUNTAIN_HEIGHT*exp(-pow(distance, 2)/(2*pow(sigma_mountain, 2)));
		}
		// reall orography can only be different from zero on land point
		if (ORO_ID == 2 && is_land[i] == 1)
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
	if (ORO_ID == 2)
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
	free(latitude_scalar);
	free(longitude_scalar);
	free(oro_unfiltered);
	
	printf("minimum orography: %lf m\n", oro[find_min_index(oro, NO_OF_SCALARS_H)]);
	printf("maximum orography: %lf m\n", oro[find_max_index(oro, NO_OF_SCALARS_H)]);
	
	// surface properties other than orography
   	double *roughness_length = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_albedo = calloc(NO_OF_SCALARS_H, sizeof(double));
	double *sfc_rho_c = calloc(NO_OF_SCALARS_H, sizeof(double));
	double c_p_water = 4184.0;
	double c_p_soil = 830.0;
	double albedo_water = 0.06;
	// setting the land surface albedo to 0.12 (compare Zdunkowski, Trautmann & Bott:
	// Radiation in the Atmosphere, 2007, p. 444)
	double albedo_soil = 0.12;
	double density_water = 1024.0;
	double density_soil = 1442.0;
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// ocean
		sfc_albedo[i] = albedo_water;
		sfc_rho_c[i] = density_water*c_p_water;
		// land
		if (is_land[i] == 1)
		{
			sfc_albedo[i] = albedo_soil;
			sfc_rho_c[i] = density_soil*c_p_soil;
		}
		// roughness length
		roughness_length[i] = 0.02;
	}
	int sfc_albedo_id, sfc_rho_c_id, roughness_length_id;
	if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
	  ERR(retval);
	if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS_H, &scalar_h_dimid)))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "z_surface", NC_DOUBLE, 1, &scalar_h_dimid, &oro_id)))
	  ERR(retval);
	if ((retval = nc_put_att_text(ncid, oro_id, "units", strlen("m"), "m")))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "sfc_albedo", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_albedo_id)))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "sfc_rho_c", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_rho_c_id)))
	  ERR(retval);
	if ((retval = nc_put_att_text(ncid, sfc_rho_c_id, "units", strlen("J/(K*m**3)"), "J/(K*m**3)")))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "is_land", NC_INT, 1, &scalar_h_dimid, &is_land_id)))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "roughness_length", NC_DOUBLE, 1, &scalar_h_dimid, &roughness_length_id)))
	  ERR(retval);
	if ((retval = nc_put_att_text(ncid, roughness_length_id, "units", strlen("m"), "m")))
	  ERR(retval);
	if ((retval = nc_enddef(ncid)))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, oro_id, &oro[0])))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, sfc_albedo_id, &sfc_albedo[0])))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, sfc_rho_c_id, &sfc_rho_c[0])))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, roughness_length_id, &roughness_length[0])))
	  ERR(retval);
	if ((retval = nc_put_var_int(ncid, is_land_id, &is_land[0])))
	  ERR(retval);
	if ((retval = nc_close(ncid)))
	  ERR(retval);
	free(roughness_length);
	free(sfc_albedo);
	free(sfc_rho_c);
	free(is_land);
	free(oro);
	
	return 0;
}











