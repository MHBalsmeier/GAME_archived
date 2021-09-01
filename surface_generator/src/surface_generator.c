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
	int OUTPUT_FILE_LENGTH = 100;
	char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
	int ORO_ID;
   	ORO_ID = strtod(argv[1], NULL);
   	if (ORO_ID < 1 || ORO_ID > 2)
   	{
   		printf("Error: oro_id must not be smaller than one or larger than 2.\n");
   		exit(1);
	}
   	double oro_rescale_factor = strtof(argv[2], NULL);
	sprintf(OUTPUT_FILE_PRE, "surface_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
	OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
	free(OUTPUT_FILE_PRE);
	char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
	sprintf(OUTPUT_FILE, "surface_files/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
	int ncid, scalar_h_dimid, oro_id, latitude_scalar_id, longitude_scalar_id;
	double *oro = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_albedo = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_c_v = malloc(NO_OF_SCALARS_H*sizeof(double));
    int GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/grids/B%dL26T41152_O0_OL23_SCVT.nc", RES_ID);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/grids/B%dL26T41152_O0_OL23_SCVT.nc", RES_ID);
	int retval;
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double distance;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        ERR(retval);
    free(GEO_PROP_FILE);
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
		int INPUT_FILE_LENGTH = 100;
		char *INPUT_FILE_PRE = malloc((INPUT_FILE_LENGTH + 1)*sizeof(char));
		sprintf(INPUT_FILE_PRE, "real/hgt.sfc.nc");
		INPUT_FILE_LENGTH = strlen(INPUT_FILE_PRE);
		free(INPUT_FILE_PRE);
		char *INPUT_FILE = malloc((INPUT_FILE_LENGTH + 1)*sizeof(char));
		sprintf(INPUT_FILE, "real/hgt.sfc.nc");
		if ((retval = nc_open(INPUT_FILE, NC_NOWRITE, &ncid)))
		    ERR(retval);
		free(INPUT_FILE);
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
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
    	// ocean
    	sfc_albedo[i] = 0.06;
    	sfc_c_v[i] = 4184;
		// setting the land surface albedo to 0.12 (compare Zdunkowski,Trautmann & Bott:
		// Radiation in the Atmosphere,2007,p. 444)
    	if (oro[i] > 5)
    	{
    		sfc_albedo[i] = 0.12;
    		sfc_c_v[i] = 0.5*4184;
    	}
    }
    int sfc_albedo_id, sfc_c_v_id;
	if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
	  ERR(retval);
	free(OUTPUT_FILE);
	if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS_H, &scalar_h_dimid)))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "z_surface", NC_DOUBLE, 1, &scalar_h_dimid, &oro_id)))
	  ERR(retval);
	if ((retval = nc_put_att_text(ncid, oro_id, "units", strlen("m"), "m")))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "sfc_albedo", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_albedo_id)))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "sfc_c_v", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_c_v_id)))
	  ERR(retval);
	if ((retval = nc_put_att_text(ncid, sfc_c_v_id, "units", strlen("J/(kg*K)"), "J/(kg*K)")))
	  ERR(retval);
	if ((retval = nc_enddef(ncid)))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, oro_id, &oro[0])))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, sfc_albedo_id, &sfc_albedo[0])))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, sfc_c_v_id, &sfc_c_v[0])))
	  ERR(retval);
	if ((retval = nc_close(ncid)))
	  ERR(retval);
	free(oro);
	free(sfc_albedo);
	free(sfc_c_v);
	return 0;
}











