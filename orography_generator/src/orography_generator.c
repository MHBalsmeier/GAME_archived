/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
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
#define UNIT "Z_SURFACE"
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define OMEGA (7.292115e-5)
#define P_0 100000.0

const double U_0 = 35;
const double ETA_0 = 0.252;
const double ETA_T = 0.2;
double T_0 = 288;
const double G = 9.80616;
double GAMMA = 0.005;
const double DELTA_T = 4.8e5;

int find_z_from_p(double, double, double *);

int main(int argc, char *argv[])
{	
	int OUTPUT_FILE_LENGTH = 100;
	char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
	int ORO_ID;
   	ORO_ID = strtod(argv[1], NULL);
   	if (ORO_ID < 1 || ORO_ID > 3)
   	{
   		printf("Error: oro_id must not be smaller than one or larger than 3.\n");
   		exit(1);
	}
	sprintf(OUTPUT_FILE_PRE, "orographies/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
	OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
	free(OUTPUT_FILE_PRE);
	char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
	sprintf(OUTPUT_FILE, "orographies/B%d_O%d_SCVT.nc", RES_ID, ORO_ID);
	int ncid, scalar_h_dimid, var_dimid, oro_id, latitude_scalar_id, longitude_scalar_id;
	double *oro, latitude;
	oro = malloc(NO_OF_SCALARS_H*sizeof(double));
    int GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/grids/B%dL26T30000_O0_OL17_SCVT.nc", RES_ID);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/grids/B%dL26T30000_O0_OL17_SCVT.nc", RES_ID);
	int scalar_index, retval;
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
	if (ORO_ID == 3)
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
    int i;
	#pragma omp parallel for private (lat_index, lon_index)
    for (i = 0; i < no_of_lat_points*no_of_lon_points; ++i)
    {
    	lat_index = i/no_of_lon_points;
    	lon_index = i - lat_index*no_of_lon_points;
    	z_in_vector[i] = z_input[lat_index][lon_index];
    }
    double weights_sum;
    int j;
	#pragma omp parallel for private(j, distance, latitude, lat_index, lon_index, weights_sum)
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		double distance_vector[no_of_lat_points*no_of_lon_points];
		int min_indices_vector[4];
		double weights_vector[4];
		if (ORO_ID == 1)
		{
            distance = calculate_distance_h(latitude_scalar[i], longitude_scalar[i], 0, 0, RADIUS);
			oro[i] = 10e3*exp(-pow(distance/100e3, 2));
		}
		if (ORO_ID == 2)
		{
			latitude = latitude_scalar[i];
			find_z_from_p(latitude, P_0, &oro[i]);
		}
		if (ORO_ID == 3)
		{
			for (j = 0; j < no_of_lat_points*no_of_lon_points; ++j)
			{
				lat_index = j/no_of_lon_points;
				lon_index = j - lat_index*no_of_lon_points;
				distance_vector[j] = calculate_distance_h(deg2rad(latitude_input[lat_index]), deg2rad(longitude_input[lon_index]), latitude_scalar[i], longitude_scalar[i], 1);
			}
			for (j = 0; j < 4; ++j)
			{
				min_indices_vector[j] = -1;
			}
			weights_sum = 0;
			for (j = 0; j < 4; ++j)
			{
				min_indices_vector[j] = find_min_index_exclude(distance_vector, no_of_lat_points*no_of_lon_points, min_indices_vector, 4);
				weights_vector[j] = 1/(distance_vector[min_indices_vector[j]] + 0.01);
				weights_sum += weights_vector[j];
			}
			oro[i] = 0;
			for (j = 0; j < 4; ++j)
			{
				oro[i] += z_in_vector[min_indices_vector[j]]*weights_vector[j]/weights_sum;
			}
			if (oro[i] < -600 || oro[i] > 5700)
				printf("Warning: value out of usual range.\n");		
		}
	}
	free(z_in_vector);
	free(z_input);
	free(latitude_input);
	free(longitude_input);
	free(latitude_scalar);
	free(longitude_scalar);
	if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
	  ERR(retval);
	free(OUTPUT_FILE);
	if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS_H, &scalar_h_dimid)))
	  ERR(retval);
	if ((retval = nc_def_var(ncid, "z_surface", NC_DOUBLE, 1, &scalar_h_dimid, &oro_id)))
	  ERR(retval);
	if ((retval = nc_put_att_text(ncid, oro_id, "units", strlen(UNIT), UNIT)))
	  ERR(retval);
	if ((retval = nc_enddef(ncid)))
	  ERR(retval);
	if ((retval = nc_put_var_double(ncid, oro_id, &oro[0])))
	  ERR(retval);
	if ((retval = nc_close(ncid)))
	  ERR(retval);
	free(oro);
	return 0;
}

int find_z_from_p(double lat, double p, double *result)
{
    double z;
    double eta = p/P_0;
    double phi, phi_bg, phi_perturb;
    double eta_v = (eta - ETA_0)*M_PI/2;
    if (eta >= ETA_T)
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G));
    else
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G)) - R_D*DELTA_T*((log(eta/ETA_T) + 137.0/60.0)*pow(ETA_T, 5) - 5*eta*pow(ETA_T, 4) + 5*pow(ETA_T, 3)*pow(eta, 2) - 10.0/3.0*pow(ETA_T, 2)*pow(eta, 3) + 5.0/4.0*ETA_T*pow(eta, 4) - 1.0/5.0*pow(eta, 5));
    phi_perturb = U_0*pow(cos(eta_v), 1.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3.0) + 10.0/63.0)*U_0*pow(cos(eta_v), 1.5) + RADIUS*OMEGA*(8.0/5.0*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3.0) - M_PI/4.0));
    phi = phi_bg + phi_perturb;
    z = phi/G;
    *result = z;
    return 0;
}














