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
#include "geos95.h"
#include "enum.h"
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
   	double oro_rescale_factor = strtof(argv[2], NULL);
   	char OUTPUT_FILE_PRE[200];
	sprintf(OUTPUT_FILE_PRE, "surface_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
   	char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
	strcpy(OUTPUT_FILE, OUTPUT_FILE_PRE);
	int ncid, scalar_h_dimid, oro_id, latitude_scalar_id, longitude_scalar_id;
	double *oro = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_albedo = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_c_v = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_rho = malloc(NO_OF_SCALARS_H*sizeof(double));
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
   	int lat_in_id, lon_in_id, z_in_id;
   	int no_of_lat_points = 73;
   	int no_of_lon_points = 144;
    float *latitude_input = malloc(no_of_lat_points*sizeof(double));
   	float *longitude_input = malloc(no_of_lon_points*sizeof(double));
    float (*z_input)[no_of_lon_points] = malloc(sizeof(float[no_of_lat_points][no_of_lon_points]));
	if (ORO_ID == 2)
	{
		if ((retval = nc_open("real/hgt.sfc.nc", NC_NOWRITE, &ncid)))
		    ERR(retval);
		if ((retval = nc_inq_varid(ncid, "lat", &lat_in_id)))
		    ERR(retval);
		if ((retval = nc_inq_varid(ncid, "lon", &lon_in_id)))
		    ERR(retval);
		if ((retval = nc_inq_varid(ncid, "hgt", &z_in_id)))
		    ERR(retval);
		if ((retval = nc_get_var_float(ncid, lat_in_id, &latitude_input[0])))
		    ERR(retval);
		if ((retval = nc_get_var_float(ncid, lon_in_id, &longitude_input[0])))
		    ERR(retval);
		if ((retval = nc_get_var_float(ncid, z_in_id, &z_input[0][0])))
		    ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
	}
    double *z_in_vector = malloc(no_of_lat_points*no_of_lon_points*sizeof(double));
    int lat_index, lon_index;
	#pragma omp parallel for private (lat_index, lon_index)
    for (int i = 0; i < no_of_lat_points*no_of_lon_points; ++i)
    {
    	lat_index = i/no_of_lon_points;
    	lon_index = i - lat_index*no_of_lon_points;
    	z_in_vector[i] = z_input[lat_index][lon_index];
    }
    // finding the number of points used for the average
    int no_of_avg_points;
    // this is the average number of points of the input dataset per model grid cell 
    double area_ratio = sqrt((no_of_lat_points*no_of_lon_points + 0.0)/NO_OF_SCALARS_H);
    int no_of_cells_for_gliding_avrg = 9;
    no_of_avg_points = (int) no_of_cells_for_gliding_avrg*area_ratio;
    if (no_of_avg_points >= 1)
    {
    	printf("number of points used for averaging: %d\n", no_of_avg_points);
    }
    else
    {
    	printf("Error: number of points used for averaging is zero.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    // executing the actual interpolation
    double weights_sum, sigma_mountain;
	#pragma omp parallel for private(distance, lat_index, lon_index, weights_sum)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		double distance_vector[no_of_lat_points*no_of_lon_points];
		int min_indices_vector[no_of_avg_points];
		double weights_vector[no_of_avg_points];
		if (ORO_ID == 1)
		{
			sigma_mountain = MOUNTAIN_FWHM/pow(8*log(2), 0.5);
            distance = calculate_distance_h(latitude_scalar[i], longitude_scalar[i], 0, 0, RADIUS);
			oro[i] = MOUNTAIN_HEIGHT*exp(-pow(distance, 2)/(2*pow(sigma_mountain, 2)));
		}
		if (ORO_ID == 2)
		{
			for (int j = 0; j < no_of_lat_points*no_of_lon_points; ++j)
			{
				lat_index = j/no_of_lon_points;
				lon_index = j - lat_index*no_of_lon_points;
				distance_vector[j] = calculate_distance_h(deg2rad(latitude_input[lat_index]), deg2rad(longitude_input[lon_index]), latitude_scalar[i], longitude_scalar[i], 1);
			}
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				min_indices_vector[j] = -1;
			}
			weights_sum = 0;
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				min_indices_vector[j] = find_min_index_exclude(distance_vector, no_of_lat_points*no_of_lon_points, min_indices_vector, no_of_avg_points);
				weights_vector[j] = 1/(pow(distance_vector[min_indices_vector[j]], 2 + EPSILON_SECURITY) + EPSILON_SECURITY);
				weights_sum += weights_vector[j];
			}
			oro[i] = 0;
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				oro[i] += oro_rescale_factor*z_in_vector[min_indices_vector[j]]*weights_vector[j]/weights_sum;
			}
			if (oro[i] < -600 || oro[i] > 5700)
			{
				printf("Warning: value out of usual range.\n");		
			}
		}
	}
	printf("minimum orography: %lf m\n", oro[find_min_index(oro, NO_OF_SCALARS_H)]);
	printf("maximum orography: %lf m\n", oro[find_max_index(oro, NO_OF_SCALARS_H)]);
	free(z_in_vector);
	free(z_input);
	free(latitude_input);
	free(longitude_input);
	free(latitude_scalar);
	free(longitude_scalar);
	
	// other surface properties
	// reading the land mask
	if (ORO_ID == 2)
	{
		char IS_LAND_FILE_PRE[200];
		sprintf(IS_LAND_FILE_PRE, "real/B%d_is_land.nc", RES_ID);
		char IS_LAND_FILE[strlen(IS_LAND_FILE_PRE) + 1];
		strcpy(IS_LAND_FILE, IS_LAND_FILE_PRE);
		int *is_land = malloc(NO_OF_SCALARS_H*sizeof(int));
		int is_land_id;
		if ((retval = nc_open(IS_LAND_FILE, NC_NOWRITE, &ncid)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "is_land", &is_land_id)))
			ERR(retval);
		if ((retval = nc_get_var_int(ncid, is_land_id, &is_land[0])))
			ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
		double c_v_water = 4184.0;
		double c_v_soil = 830.0;
		double albedo_water = 0.06;
		double albedo_soil = 0.12;
		double density_water = 1024.0;
		double density_soil = 1442.0;
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			// ocean
			sfc_albedo[i] = albedo_water;
			sfc_c_v[i] = c_v_water;
			sfc_rho[i] = density_water;
			if (is_land[i] == 1)
			{
				// setting the land surface albedo to 0.12 (compare Zdunkowski,Trautmann & Bott:
				// Radiation in the Atmosphere,2007,p. 444)
				sfc_albedo[i] = albedo_soil;
				sfc_c_v[i] = c_v_soil;
				sfc_rho[i] = density_soil;
			}
		}
		int sfc_albedo_id, sfc_c_v_id, sfc_density_id;
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
		if ((retval = nc_def_var(ncid, "sfc_density", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_density_id)))
		  ERR(retval);
		if ((retval = nc_put_att_text(ncid, sfc_density_id, "units", strlen("kg/(m**3)"), "kg/(m**3)")))
		  ERR(retval);
		if ((retval = nc_def_var(ncid, "sfc_c_v", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_c_v_id)))
		  ERR(retval);
		if ((retval = nc_put_att_text(ncid, sfc_c_v_id, "units", strlen("J/(kg*K)"), "J/(kg*K)")))
		  ERR(retval);
		if ((retval = nc_def_var(ncid, "is_land", NC_INT, 1, &scalar_h_dimid, &is_land_id)))
		  ERR(retval);
		if ((retval = nc_enddef(ncid)))
		  ERR(retval);
		if ((retval = nc_put_var_double(ncid, oro_id, &oro[0])))
		  ERR(retval);
		if ((retval = nc_put_var_double(ncid, sfc_albedo_id, &sfc_albedo[0])))
		  ERR(retval);
		if ((retval = nc_put_var_double(ncid, sfc_density_id, &sfc_rho[0])))
		  ERR(retval);
		if ((retval = nc_put_var_double(ncid, sfc_c_v_id, &sfc_c_v[0])))
		  ERR(retval);
		if ((retval = nc_put_var_int(ncid, is_land_id, &is_land[0])))
		  ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
		free(oro);
		free(sfc_albedo);
		free(sfc_rho);
		free(sfc_c_v);
		free(is_land);
    }
	return 0;
}











