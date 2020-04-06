#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "/usr/src/eccodes/include/eccodes.h"
#include "/lib/geos/include/geos.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0
#define OMEGA (7.292115e-5)
#define C_P 1005.0

const double TROPO_HEIGHT = 12e3;
const double T_SFC = 273.15 + 15;
double T_0 = 288;
const double TEMP_GRADIENT = -0.65/100;
double GAMMA = 0.005;
const double G = 9.80616;
const double DELTA_T = 4.8e5;
const double ETA_T = 0.2;
const double U_0 = 35;
const double ETA_0 = 0.252;
const double TOA = 30000;
const short MODE = 2;
const short ORO_ID = 0;

enum grid_integers {
RES_ID = 4,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = 6,
NUMBER_OF_LEVELS = NUMBER_OF_LAYERS + 1,
NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_VECTORS_V = NUMBER_OF_SCALARS_H,
NUMBER_OF_H_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H,
NUMBER_OF_V_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_VECTORS_V,
NUMBER_OF_VECTORS_PER_LAYER = NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V,
NUMBER_OF_TRIANGLES = (int) (NUMBER_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NUMBER_OF_SCALARS = NUMBER_OF_SCALARS_H*NUMBER_OF_LAYERS,
NUMBER_OF_VECTORS = NUMBER_OF_H_VECTORS + NUMBER_OF_V_VECTORS,
NUMBER_OF_DUAL_SCALARS_H = NUMBER_OF_TRIANGLES,
NUMBER_OF_DUAL_VECTORS_H = 3*NUMBER_OF_TRIANGLES/2,
NUMBER_OF_DUAL_VECTORS_V = NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_H_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_VECTORS_H,
NUMBER_OF_DUAL_V_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_V,
NUMBER_OF_DUAL_VECTORS_PER_LAYER = NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_VECTORS_V,
NUMBER_OF_DUAL_SCALARS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_VECTORS = NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_DUAL_V_VECTORS,
TRIANGLES_PER_FACE = NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

int find_pressure_value(double, double, double *);
int find_z_from_p(double, double, double *);

int main(int argc, char *argv[])
{
    int test_id = 3;
    double *direction = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *latitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *latitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *z_scalar = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NUMBER_OF_VECTORS*sizeof(double));
    int ncid, retval;
    short GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/nc_files/B%dL%dT%d_M%d_O%d.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, MODE, ORO_ID);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/nc_files/B%dL%dT%d_M%d_O%d.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, MODE, ORO_ID);
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    free(GEO_PROP_FILE);
    int direction_id, latitude_scalar_id, longitude_scalar_id, latitude_vector_id, longitude_vector_id, z_scalar_id, z_vector_id;
    if ((retval = nc_inq_varid(ncid, "direction", &direction_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_scalar", &z_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector", &z_vector_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &direction[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &latitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &longitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &z_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &z_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    const int START_SECTION = 4;
    char *SAMPLE_FILE_SCALAR = "grib_files/scalar_field_blueprint_res_id_2.grb2";
    char *SAMPLE_FILE_VECTOR = "grib_files/vector_field_blueprint_res_id_2.grb2";
    FILE *SAMPLE_SCALAR;
    FILE *SAMPLE_VECTOR;
    int err = 0;
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    short OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "grib_files/test_%d_B%dL%dT%d_M%d_O%d.grb2", test_id, RES_ID, NUMBER_OF_LAYERS, (int) TOA, MODE, ORO_ID);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "grib_files/test_%d_B%dL%dT%d_M%d_O%d.grb2", test_id, RES_ID, NUMBER_OF_LAYERS, (int) TOA, MODE, ORO_ID);
    codes_handle *handle_pot_temperature = NULL;
    codes_handle *handle_density = NULL;
    codes_handle *handle_wind_h = NULL;
    codes_handle *handle_wind_v = NULL;
    double *pressure = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *pot_temperature = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *temperature = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *rho = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *wind_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_v = malloc(NUMBER_OF_VECTORS_V*sizeof(double));
    const double TROPO_TEMP = T_SFC + TROPO_HEIGHT*TEMP_GRADIENT;
    int scalar_index;
    double z_height;
    handle_pot_temperature = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_density = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_VECTOR = fopen(SAMPLE_FILE_VECTOR, "r");
    handle_wind_h = codes_handle_new_from_file(NULL, SAMPLE_VECTOR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_VECTOR);
    double lat, lon;
    double u, eta, eta_v, T_perturb, distance, pressure_value;
    double u_p = 1.0;
    double distance_scale = SEMIMAJOR/10;
    double lat_perturb = 2*M_PI/9;
    double lon_perturb = M_PI/9;
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            lat = latitude_scalar[j];
            lon = longitude_scalar[j];
            z_height = z_scalar[i*NUMBER_OF_SCALARS_H + j];
            if (test_id == 0 || test_id == 1)
            {
                if (z_height < TROPO_HEIGHT)
                {
                    temperature[j] = T_SFC + z_height*TEMP_GRADIENT;
                    pressure[j] = P_0*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G/(R_D*TEMP_GRADIENT));
                }
                else
                {
                    temperature[j] = TROPO_TEMP;
                    pressure[j] = P_0*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height - TROPO_HEIGHT)/(R_D*TROPO_TEMP));
                }
            }
            if (test_id == 2 || test_id == 3)
            {
                find_pressure_value(lat, z_height, &pressure_value);
                pressure[j] = pressure_value;
                eta = pressure[j]/P_0;
                eta_v = (eta - ETA_0)*M_PI/2;
                T_perturb = 3.0/4.0*eta*M_PI*U_0/R_D*sin(eta_v)*pow(cos(eta_v), 0.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3.0) + 10.0/63.0)*2*U_0*pow(cos(eta_v), 1.5) + SEMIMAJOR*OMEGA*(8.0/5.0*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3.0) - M_PI/4.0));
                if (eta >= ETA_T)
                    temperature[j] = T_0*pow(eta, R_D*GAMMA/G) + T_perturb;
                else
                    temperature[j] = T_0*pow(eta, R_D*GAMMA/G) + DELTA_T*pow(ETA_T - eta, 5) + T_perturb;
            }
            rho[j] = pressure[j]/(R_D*temperature[j]);
            pot_temperature[j] = temperature[j]*pow(P_0/pressure[j], R_D/C_P);
        }
        if (retval = codes_set_long(handle_pot_temperature, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "parameterCategory", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "parameterNumber", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_pot_temperature, "values", pot_temperature, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        if (i == 0)
            codes_write_message(handle_pot_temperature, OUTPUT_FILE, "w");
        else
            codes_write_message(handle_pot_temperature, OUTPUT_FILE, "a");
        if (retval = codes_set_long(handle_density, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "parameterCategory", 3))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "parameterNumber", 10))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_density, "values", rho, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        codes_write_message(handle_density, OUTPUT_FILE, "a");
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = z_vector[NUMBER_OF_VECTORS_V + j + i*NUMBER_OF_VECTORS_PER_LAYER];
            if (test_id == 0)
                wind_h[j] = 0;
            if (test_id == 1)
            {
                if (i == NUMBER_OF_LAYERS - 2 && j == NUMBER_OF_SCALARS_H/2)
                    wind_h[j] = 10;
                else
                    wind_h[j] = 0;
            }
            if (test_id == 2 || test_id == 3)
            {
                find_pressure_value(lat, z_height, &pressure_value);
                eta = pressure_value/P_0;
                eta_v = (eta - ETA_0)*M_PI/2; 
                u = U_0*pow(cos(eta_v), 1.5)*pow(sin(2*lat), 2);
                if (test_id == 3)
                {
                    distance = calculate_distance_h(lat, lon, lat_perturb, lon_perturb, SEMIMAJOR);
                    u += u_p*exp(-pow(distance/distance_scale, 2));
                }
                wind_h[j] = u*cos(direction[j]);
            }
        }
        if (retval = codes_set_long(handle_wind_h, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "parameterCategory", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "parameterNumber", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_wind_h, "values", wind_h, NUMBER_OF_VECTORS_H))
            ECCERR(retval);
        codes_write_message(handle_wind_h, OUTPUT_FILE, "a");
    }
    free(z_scalar);
    codes_handle_delete(handle_pot_temperature);
    codes_handle_delete(handle_density);
    codes_handle_delete(handle_wind_h);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_wind_v = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    for (int i = 0; i < NUMBER_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS_V; ++j)
        {
            lat = latitude_scalar[j];
            lon = longitude_scalar[j];
            z_height = z_vector[j + i*NUMBER_OF_VECTORS_PER_LAYER];
            if (test_id == 0 || test_id == 1 || test_id == 2 || test_id == 3)
                wind_v[j] = 0;
        }
        if (retval = codes_set_long(handle_wind_v, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "parameterCategory", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "parameterNumber", 9))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_wind_v, "values", wind_v, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        codes_write_message(handle_wind_v, OUTPUT_FILE, "a");
    }
    free(z_vector);
    free(latitude_vector);
    free(longitude_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    free(direction);
    free(OUTPUT_FILE);
    free(pressure);
    free(pot_temperature);
    free(temperature);
    free(rho);
    free(wind_h);
    free(wind_v);
    return 0;
}

int find_pressure_value(double lat, double z_height, double *result)
{
    double p = P_0/2;
    double precision = 0.0001;
    double z;
    double current_max = P_0;
    double current_min = 0;
    find_z_from_p(lat, p, &z);
    while (fabs(z - z_height) > precision)
    {
        if (z < z_height)
        {
            current_max = p;
            p = 0.5*(current_min + p);
        }
        else
        {
            current_min = p;
            p = 0.5*(current_max + p);
        }
        find_z_from_p(lat, p, &z);
    }
    *result = p;
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
    phi_perturb = U_0*pow(cos(eta_v), 1.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3.0) + 10.0/63.0)*U_0*pow(cos(eta_v), 1.5) + SEMIMAJOR*OMEGA*(8.0/5.0*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3.0) - M_PI/4.0));
    phi = phi_bg + phi_perturb;
    z = phi/G;
    *result = z;
    return 0;
}

















