/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include "geos95.h"
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

const int MODE = 2;
const int ORO_ID = 2;
const double U_0 = 35;
const double ETA_0 = 0.252;
const double ETA_T = 0.2;
double T_0 = 288;
const double G = 9.80616;
double GAMMA = 0.005;
const double DELTA_T = 4.8e5;

/*
ORO_ID:
0	sphere
1	sphere with Gaussian mountain at 0 / 0, H = 10 km
2	JW test
*/

int find_z_from_p(double, double, double *);

enum grid_integers {
RES_ID = 4,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = 6,
NUMBER_OF_ORO_LAYERS = 4,
NUMBER_OF_LEVELS = NUMBER_OF_LAYERS + 1,
NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_VECTORS_V = NUMBER_OF_SCALARS_H,
NUMBER_OF_VECTORS_PER_LAYER = NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V,
NUMBER_OF_TRIANGLES = (int) (NUMBER_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NUMBER_OF_DUAL_SCALARS_H = NUMBER_OF_TRIANGLES,
NUMBER_OF_DUAL_VECTORS_H = 3*NUMBER_OF_TRIANGLES/2,
NUMBER_OF_DUAL_VECTORS_V = NUMBER_OF_DUAL_SCALARS_H,
TRIANGLES_PER_FACE = NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

int main(int argc, char *argv[])
{
	int OUTPUT_FILE_LENGTH = 100;
	char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
	sprintf(OUTPUT_FILE_PRE, "nc_files/B%d_M%d_O%d.nc", RES_ID, MODE, ORO_ID);
	OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
	free(OUTPUT_FILE_PRE);
	char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
	sprintf(OUTPUT_FILE, "nc_files/B%d_M%d_O%d.nc", RES_ID, MODE, ORO_ID);
	int ncid, scalar_h_dimid, var_dimid, oro_id, latitude_scalar_id, longitude_scalar_id;
	double *oro, latitude;
	oro = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    int GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/nc_files/B%dL26T30000_M%d_O0_OL17.nc", RES_ID, MODE);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/nc_files/B%dL26T30000_M%d_O0_OL17.nc", RES_ID, MODE);
	int scalar_index, retval;
    double *latitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
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
	for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
	{	
		if (ORO_ID == 0)
			oro[i] = 0;
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
	}
	free(latitude_scalar);
	free(longitude_scalar);
	if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
	  ERR(retval);
	free(OUTPUT_FILE);
	if ((retval = nc_def_dim(ncid, "scalar_index", NUMBER_OF_SCALARS_H, &scalar_h_dimid)))
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














