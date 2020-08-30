/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
The grid generation procedure is manged from this file. Memory allocation and IO is done here, for the rest, functions are called residing in individual files.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include "grid_generator.h"
#include "enum.h"
#include "geos95.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define OMEGA (7.292115e-5)
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define C_P 1005.0
#define P_0 100000.0
#define GRAVITY_MEAN_SFC_ABS 9.80616
#define RESET "\033[0m"
#define BLACK "\033[30m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"
#define BOLDBLACK "\033[1m\033[30m"
#define BOLDRED "\033[1m\033[31m"
#define BOLDGREEN "\033[1m\033[32m"
#define BOLDYELLOW "\033[1m\033[33m"
#define BOLDBLUE "\033[1m\033[34m"
#define BOLDMAGENTA "\033[1m\033[35m"
#define BOLDCYAN "\033[1m\033[36m"
#define BOLDWHITE "\033[1m\033[37m"

const double TOA = 30000.0;
const double ORTH_CRITERION_DEG = 89.99;

int main(int argc, char *argv[])
{
    int ORO_ID;
   	ORO_ID = strtod(argv[1], NULL);
    int OPTIMIZE_BOOL;
   	OPTIMIZE_BOOL = strtod(argv[2], NULL);
    int N_ITERATIONS;
   	N_ITERATIONS = strtod(argv[3], NULL);
    int USE_SCALAR_H_FILE;
   	USE_SCALAR_H_FILE = strtod(argv[4], NULL);
    int len = strlen(argv[5]);
    char *SCALAR_H_FILE = malloc((len + 1)*sizeof(char));
    strcpy(SCALAR_H_FILE, argv[5]);
	if (NO_OF_ORO_LAYERS >= NO_OF_LAYERS)
	{
		printf("It is NO_OF_ORO_LAYERS >= NO_OF_LAYERS.\n");
		exit(1);
	}
    char OUTPUT_FILE_PRE[200];
    char STATISTICS_FILE_PRE[200];
    if (OPTIMIZE_BOOL == 1)
    {
    	sprintf(OUTPUT_FILE_PRE, "grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    	sprintf(STATISTICS_FILE_PRE, "statistics/B%dL%dT%d_O%d_OL%d_SCVT.txt", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
	}
	else
	{
		if (USE_SCALAR_H_FILE == 1)
		{
			if (SCALAR_H_FILE[strlen(SCALAR_H_FILE) - 1 - 6] == 'S' && SCALAR_H_FILE[strlen(SCALAR_H_FILE) - 1 - 5] == 'C' && 
			SCALAR_H_FILE[strlen(SCALAR_H_FILE) - 1 - 4] == 'V' && SCALAR_H_FILE[strlen(SCALAR_H_FILE) - 1 - 3] == 'T')
			{
    			sprintf(OUTPUT_FILE_PRE, "grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    			sprintf(STATISTICS_FILE_PRE, "statistics/B%dL%dT%d_O%d_OL%d_SCVT.txt", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
			}
    		else
    		{
    			sprintf(OUTPUT_FILE_PRE, "grids/B%dL%dT%d_O%d_OL%d.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    			sprintf(STATISTICS_FILE_PRE, "statistics/B%dL%dT%d_O%d_OL%d.txt", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
			}
		}
		else
		{
    		sprintf(OUTPUT_FILE_PRE, "grids/B%dL%dT%d_O%d_OL%d.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    		sprintf(STATISTICS_FILE_PRE, "statistics/B%dL%dT%d_O%d_OL%d.txt", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
		}
	}
    char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
    char STATISTICS_FILE[strlen(STATISTICS_FILE_PRE) + 1];
    strcpy(OUTPUT_FILE, OUTPUT_FILE_PRE);
    strcpy(STATISTICS_FILE, STATISTICS_FILE_PRE);
	printf("Output will be written to file %s.\n", OUTPUT_FILE);
    double *latitude_ico = malloc(12*sizeof(double));
    double *longitude_ico = malloc(12*sizeof(double));
    int edge_vertices[NO_OF_EDGES][2];
    int face_vertices[20][3];
    int face_edges[20][3];
    int face_edges_reverse[20][3];
    printf("Building icosahedron ... ");
	build_icosahedron(latitude_ico, longitude_ico, edge_vertices, face_vertices, face_edges, face_edges_reverse);
    printf(GREEN "finished.\n" RESET);
    printf("Allocating memory ... ");
    double *x_unity = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *y_unity = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_unity = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_scalar = malloc(NO_OF_SCALARS*sizeof(double));
    double *gravity_potential = malloc(NO_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NO_OF_VECTORS*sizeof(double));
    double *normal_distance = malloc(NO_OF_VECTORS*sizeof(double));
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *direction = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *volume = malloc(NO_OF_SCALARS*sizeof(double));
    double *area = malloc(NO_OF_VECTORS*sizeof(double));
    double *trsk_modified_weights = calloc(10*NO_OF_VECTORS_H, sizeof(double));
    double *recov_ver_weight = malloc(6*NO_OF_LEVELS*NO_OF_SCALARS_H*sizeof(double));
    double *latitude_scalar_dual = malloc(NO_OF_DUAL_SCALARS_H*sizeof(double));
    double *longitude_scalar_dual = malloc(NO_OF_DUAL_SCALARS_H*sizeof(double));
    double *z_scalar_dual = malloc(NO_OF_DUAL_SCALARS*sizeof(double));
    double *z_vector_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *normal_distance_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *direction_dual = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *area_dual_pre = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *area_dual = malloc((NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS)*sizeof(double));
    double *f_vec = malloc(3*NO_OF_VECTORS_H*sizeof(double));
    double *triangle_face_unit_sphere = malloc(NO_OF_DUAL_SCALARS_H*sizeof(double));
    double *pent_hex_face_unity_sphere = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *rel_on_line_dual = malloc(NO_OF_VECTORS_H*sizeof(double));
	double *z_surface = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *e_kin_weights = malloc(8*NO_OF_SCALARS*sizeof(double));
	double *volume_ratios = malloc(2*NO_OF_SCALARS*sizeof(double));
	double *slope = malloc(NO_OF_VECTORS*sizeof(double));
    double *recov_primal2dual_weights = malloc(2*NO_OF_DUAL_H_VECTORS*sizeof(double));
    double *density_to_rhombus_weights = malloc(4*NO_OF_VECTORS_H*sizeof(double));
    int *to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *from_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *trsk_modified_velocity_indices = calloc(10*NO_OF_VECTORS_H, sizeof(int));
    int *trsk_modified_curl_indices = calloc(10*NO_OF_VECTORS_H, sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *vorticity_indices_pre = malloc(3*NO_OF_DUAL_SCALARS_H*sizeof(int));
    int *vorticity_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *h_curl_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *to_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *from_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_signs_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *vorticity_signs_pre = malloc(3*NO_OF_DUAL_SCALARS_H*sizeof(int));
    int *vorticity_signs = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *h_curl_signs = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *density_to_rhombus_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    printf(GREEN "finished.\n" RESET);
    int retval, z_surface_id;
    printf("Reading orography data ... ");
	set_orography(RES_ID, ORO_ID, z_surface);
    printf(GREEN "finished.\n" RESET);
    if (NO_OF_VECTORS_H != NO_OF_VECTORS_H)
    {
        printf("It is NO_OF_VECTORS_H != NO_OF_VECTORS_H.\n");
    }
    printf("Establishing horizontal grid structure ... \n");
    if (USE_SCALAR_H_FILE == 0)
    {
    	// Here, the positions of the horizontal generators, i.e. the horizontal scalar points are determined.
    	generate_horizontal_generators(latitude_ico, longitude_ico, latitude_scalar, longitude_scalar, x_unity, y_unity, z_unity, face_edges_reverse, face_edges, face_vertices);
    	// By setting the from_index and to_index arrrays, the discrete positions of the vector points are determined.
    	set_from_to_index(from_index, to_index, face_edges, face_edges_reverse, face_vertices, edge_vertices);
    	// By setting the from_index_dual and to_index_dual arrrays, the discrete positions of the dual scalar points are determined.
		set_from_to_index_dual(from_index_dual, to_index_dual, face_edges, face_edges_reverse);
    }
    else
    {
    	read_horizontal_explicit(latitude_scalar, longitude_scalar, from_index, to_index, from_index_dual, to_index_dual, SCALAR_H_FILE);
    }
	if (OPTIMIZE_BOOL == 1)
	{
		optimize_to_scvt(latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, N_ITERATIONS, face_edges, face_edges_reverse, face_vertices, edge_vertices, adjacent_vector_indices_h, from_index_dual, to_index_dual);
	}
	// Calculation of the horizontal coordinates of the dual scalar points. The dual scalar points are the vertices of the Voronoi mesh of the primal grid.
	set_scalar_h_dual_coords(latitude_scalar_dual, longitude_scalar_dual, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices, edge_vertices);
	// Calculation of the horizontal coordinates of the vector points. The vector points are in the middle between the neighbouring scalar points.
	set_vector_h_doubles(from_index, to_index, latitude_scalar, longitude_scalar, latitude_vector, longitude_vector, direction);
	// Same as before but for dual vectors. They have the same positions as the primal vectors.
	set_dual_vector_h_doubles(latitude_scalar_dual, latitude_vector, direction_dual, longitude_vector, to_index_dual, from_index_dual, longitude_scalar_dual, rel_on_line_dual);
	// Setting the Coriolis vector.
    set_f_vec(latitude_vector, direction, direction_dual, f_vec);
    printf(GREEN "Horizontal grid structure determined.\n" RESET);
    calc_triangle_face_unity(triangle_face_unit_sphere, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices, edge_vertices);
	calc_vorticity_indices_pre(from_index_dual, to_index_dual, direction, direction_dual, vorticity_indices_pre, ORTH_CRITERION_DEG, vorticity_signs_pre);
	check_for_orthogonality(direction, direction_dual, ORTH_CRITERION_DEG);
	calc_adjacent_vector_indices_h(from_index, to_index, adjacent_signs_h, adjacent_vector_indices_h);
	calc_cell_face_unity(pent_hex_face_unity_sphere, latitude_scalar_dual, longitude_scalar_dual, adjacent_vector_indices_h, vorticity_indices_pre);
	// building the vertical grid
   	double z_oro_off = TOA*(NO_OF_ORO_LAYERS + 0.0)/NO_OF_LAYERS;
	double z_vertical_vector_pre[NO_OF_LAYERS + 1];
	determine_z_scalar(z_scalar, z_vertical_vector_pre, z_surface, z_oro_off, TOA);
	set_z_vector_and_normal_distance(z_vector, z_surface, z_scalar, normal_distance, latitude_scalar, longitude_scalar, from_index, to_index, TOA);
	printf("Mapping horizontal areas from unit sphere to model levels ... ");
	map_area_to_sphere(area, z_vector, pent_hex_face_unity_sphere);
    printf(GREEN "finished.\n" RESET);
    printf("Calculating grid box volumes ... ");
	set_volume(volume, z_scalar_dual, z_vector, area, from_index, to_index, TOA, z_surface, vorticity_indices_pre);
    printf(GREEN "finished.\n" RESET);
	printf("Determining vector z coordinates of the dual grid and distances of the dual grid ... ");
	calc_z_vector_dual_and_normal_distance_dual(z_vector_dual, normal_distance_dual, z_scalar_dual, TOA, z_surface, from_index, to_index, z_vector, from_index_dual, to_index_dual, latitude_scalar_dual, longitude_scalar_dual, vorticity_indices_pre);
    printf(GREEN "finished.\n" RESET);
    printf("Determining coordinate slopes ...");
    calc_slopes(z_scalar, from_index, to_index, normal_distance, slope);
    printf(GREEN "finished.\n" RESET);
    printf("Calculating dual areas, pre version ... ");
	calc_area_dual_pre(area_dual_pre, z_vector_dual, normal_distance, z_vector, from_index, to_index, z_surface, triangle_face_unit_sphere, TOA);
    printf(GREEN "finished.\n" RESET);
    printf("Calculating vertical faces, pre version ... ");
	calculate_vertical_faces(area, z_vector_dual, normal_distance_dual, TOA);
    printf(GREEN "finished.\n" RESET);
    printf("Modifying dual areas for rhombus vorticity calculation ... ");
	modify_area_dual(area_dual, z_vector, from_index_dual, to_index_dual, area_dual_pre, z_vector_dual);
    printf(GREEN "finished.\n" RESET);
    // more advanced stuff: tangential vector reconstruction and kinetic energy
    printf("Calculating kinetic energy indices and weights and related things ... ");
	calc_kinetic_energy_and_related(latitude_scalar, longitude_scalar, e_kin_weights, volume, to_index, from_index, area_dual_pre, area, z_scalar, z_vector, adjacent_vector_indices_h, latitude_vector, longitude_vector, latitude_scalar_dual, longitude_scalar_dual, to_index_dual, from_index_dual, z_vector_dual, volume_ratios, recov_primal2dual_weights);
    set_recov_ver(adjacent_vector_indices_h, direction, direction_dual, latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, from_index_dual, to_index_dual, pent_hex_face_unity_sphere, recov_ver_weight, ORTH_CRITERION_DEG, normal_distance, z_vector, z_vector_dual, TOA);
    printf(GREEN "finished.\n" RESET);
    // modified TRSK
    printf("Calculating Coriolis indices and weights ... ");
	calc_coriolis_weights(from_index_dual, to_index_dual, trsk_modified_curl_indices, normal_distance, normal_distance_dual, to_index, area, z_scalar, latitude_scalar, longitude_scalar, latitude_vector, longitude_vector, latitude_scalar_dual, longitude_scalar_dual, trsk_modified_weights,  trsk_modified_velocity_indices, from_index, adjacent_vector_indices_h, direction, direction_dual, ORTH_CRITERION_DEG, z_vector, z_vector_dual);
    printf(GREEN "finished.\n" RESET);
    // setting indices related to the curl of a vector field
    printf("Setting horizontal curl indices ... ");
	set_horizontal_curl_indices(direction_dual, direction, h_curl_indices, from_index, to_index, ORTH_CRITERION_DEG, h_curl_signs);
    printf(GREEN "finished.\n" RESET);
    printf("Setting vorticity indices and rhombus density interpolation indices and weights ... ");
	set_vertical_vorticity_stuff(vorticity_indices_pre, vorticity_signs_pre, from_index_dual, to_index_dual, vorticity_indices, vorticity_signs, density_to_rhombus_indices, from_index, to_index, area_dual, z_vector, latitude_scalar_dual, longitude_scalar_dual, density_to_rhombus_weights, latitude_vector, longitude_vector, latitude_scalar, longitude_scalar);
    printf(GREEN "finished.\n" RESET);
    printf("Setting gravity potential ... ");
	set_gravity_potential(z_scalar, gravity_potential, GRAVITY_MEAN_SFC_ABS);
    printf(GREEN "finished.\n" RESET);
    // A statistics file is created to compare the fundamental statistical properties of the grid with the literature.
	write_statistics_file(pent_hex_face_unity_sphere, normal_distance, normal_distance_dual, STATISTICS_FILE);
    int latitude_scalar_id, longitude_scalar_id, direction_id, latitude_vector_id, longitude_vector_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_id, z_vector_id, normal_distance_id, volume_id, area_id, trsk_modified_weights_id, recov_ver_weight_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, trsk_modified_velocity_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, f_vec_dimid, scalar_dimid, scalar_h_dimid, scalar_dual_h_dimid, vector_dimid, scalar_h_dimid_6, vector_h_dimid, vector_h_dimid_11, vector_h_dimid_10, vector_h_dimid_2, vector_h_dimid_4, vector_v_dimid_6, vector_dual_dimid, vector_dual_h_dimid, vector_dual_v_dimid_3, gravity_potential_id, scalar_dual_h_dimid_3, vector_dual_area_dimid, e_kin_weights_id, scalar_8_dimid, slope_id, scalar_2_dimid, volume_ratios_id, recov_primal2dual_weights_id, vector_h_dual_dimid_2, density_to_rhombus_indices_id, density_to_rhombus_weights_id, ncid_g_prop;
    printf("Starting to write to output file ... ");
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid_g_prop)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_8_index", 8*NO_OF_SCALARS, &scalar_8_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_2_index", 2*NO_OF_SCALARS, &scalar_2_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_index", NO_OF_SCALARS_H, &scalar_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_h_index", NO_OF_DUAL_SCALARS_H, &scalar_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_h_3_index", 3*NO_OF_DUAL_SCALARS_H, &scalar_dual_h_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_index", NO_OF_VECTORS_H, &vector_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_6_index", 6*NO_OF_SCALARS_H, &scalar_h_dimid_6)))
        ERR(retval);
	if ((retval = nc_def_dim(ncid_g_prop, "vector_h_10_index", 10*NO_OF_VECTORS_H, &vector_h_dimid_10)))
	    ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_11_index", 11*NO_OF_VECTORS_H, &vector_h_dimid_11)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_2_index", 2*NO_OF_VECTORS_H, &vector_h_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_4_index", 4*NO_OF_VECTORS_H, &vector_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_6_index", 6*NO_OF_LEVELS*NO_OF_SCALARS_H, &vector_v_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_dual", NO_OF_VECTORS_H, &vector_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "f_vec_index", 3*NO_OF_VECTORS_H, &f_vec_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual", NO_OF_DUAL_VECTORS, &vector_dual_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual_area", NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS, &vector_dual_area_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_2_dual", 2*NO_OF_DUAL_H_VECTORS, &vector_h_dual_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_v_3_index", 3*NO_OF_DUAL_SCALARS_H, &vector_dual_v_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_scalar", NC_DOUBLE, 1, &scalar_h_dimid, &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_scalar", NC_DOUBLE, 1, &scalar_h_dimid, &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_scalar_dual", NC_DOUBLE, 1, &scalar_dual_h_dimid, &latitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_scalar_dual", NC_DOUBLE, 1, &scalar_dual_h_dimid, &longitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_scalar", NC_DOUBLE, 1, &scalar_dimid, &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_scalar_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "gravity_potential", NC_DOUBLE, 1, &scalar_dimid, &gravity_potential_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, gravity_potential_id, "units", strlen("m^2/s^2"), "m^2/s^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_surface", NC_DOUBLE, 1, &scalar_h_dimid, &z_surface_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_surface_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector", NC_DOUBLE, 1, &vector_dimid, &z_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "slope", NC_DOUBLE, 1, &vector_dimid, &slope_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance", NC_DOUBLE, 1, &vector_dimid, &normal_distance_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "volume", NC_DOUBLE, 1, &scalar_dimid, &volume_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, volume_id, "units", strlen("m^3"), "m^3")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "area", NC_DOUBLE, 1, &vector_dimid, &area_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, area_id, "units", strlen("m^2"), "m^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_weights", NC_DOUBLE, 1, &vector_h_dimid_10, &trsk_modified_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector_dual", NC_DOUBLE, 1, &vector_dual_dimid, &z_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance_dual", NC_DOUBLE, 1, &vector_dual_dimid, &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "area_dual", NC_DOUBLE, 1, &vector_dual_area_dimid, &area_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, area_dual_id, "units", strlen("m^2"), "m^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "f_vec", NC_DOUBLE, 1, &f_vec_dimid, &f_vec_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, f_vec_id, "units", strlen("1/s"), "1/s")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "direction", NC_DOUBLE, 1, &vector_h_dimid, &direction_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_vector", NC_DOUBLE, 1, &vector_h_dimid, &latitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_vector", NC_DOUBLE, 1, &vector_h_dimid, &longitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "e_kin_weights", NC_DOUBLE, 1, &scalar_8_dimid, &e_kin_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "volume_ratios", NC_DOUBLE, 1, &scalar_2_dimid, &volume_ratios_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_primal2dual_weights", NC_DOUBLE, 1, &vector_h_dual_dimid_2, &recov_primal2dual_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "density_to_rhombus_weights", NC_DOUBLE, 1, &vector_h_dimid_4, &density_to_rhombus_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index", NC_INT, 1, &vector_h_dimid, &from_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index", NC_INT, 1, &vector_h_dimid, &to_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index_dual", NC_INT, 1, &vector_h_dimid, &from_index_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index_dual", NC_INT, 1, &vector_h_dimid, &to_index_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_h", NC_INT, 1, &scalar_h_dimid_6, &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices", NC_INT, 1, &vector_h_dimid_4, &vorticity_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_indices", NC_INT, 1, &vector_h_dimid_4, &h_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_velocity_indices", NC_INT, 1, &vector_h_dimid_10, &trsk_modified_velocity_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_curl_indices", NC_INT, 1, &vector_h_dimid_10, &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_signs_h", NC_INT, 1, &scalar_h_dimid_6, &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs", NC_INT, 1, &vector_h_dimid_4, &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_signs", NC_INT, 1, &vector_h_dimid_4, &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "density_to_rhombus_indices", NC_INT, 1, &vector_h_dimid_4, &density_to_rhombus_indices_id)))
        ERR(retval);
    if ((retval = nc_enddef(ncid_g_prop)))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_scalar_id, &latitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_scalar_id, &longitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_scalar_dual_id, &latitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_scalar_dual_id, &longitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_scalar_id, &z_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, gravity_potential_id, &gravity_potential[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_surface_id, &z_surface[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, volume_id, &volume[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_id, &area[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, e_kin_weights_id, &e_kin_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, volume_ratios_id, &volume_ratios[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, trsk_modified_weights_id, &trsk_modified_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_weight_id, &recov_ver_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_dual_id, &z_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, direction_id, &direction[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_vector_id, &latitude_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_vector_id, &longitude_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, slope_id, &slope[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_primal2dual_weights_id, &recov_primal2dual_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, density_to_rhombus_weights_id, &density_to_rhombus_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, from_index_id, &from_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, to_index_id, &to_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, from_index_dual_id, &from_index_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, to_index_dual_id, &to_index_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vorticity_indices_id, &vorticity_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, h_curl_indices_id, &h_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, trsk_modified_velocity_indices_id, &trsk_modified_velocity_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, trsk_modified_curl_indices_id, &trsk_modified_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, h_curl_signs_id, &h_curl_signs[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, density_to_rhombus_indices_id, &density_to_rhombus_indices[0])))
        ERR(retval);
    if ((retval = nc_close(ncid_g_prop)))
        ERR(retval);
    printf(GREEN "finished.\n" RESET);
    free(SCALAR_H_FILE);
    free(x_unity);
    free(y_unity);
    free(z_unity);
    free(pent_hex_face_unity_sphere);
    free(triangle_face_unit_sphere);
    free(direction_dual);
    free(density_to_rhombus_weights);
    free(density_to_rhombus_indices);
    free(recov_primal2dual_weights);
    free(slope);
    free(rel_on_line_dual);
    free(volume_ratios);
    free(e_kin_weights);
    free(gravity_potential);
    free(latitude_vector);
    free(longitude_vector);
    free(direction);
    free(latitude_scalar);
    free(longitude_scalar);
    free(z_scalar);
    free(z_vector);
    free(normal_distance);
    free(volume);
    free(area);
    free(trsk_modified_weights);
    free(recov_ver_weight);
    free(latitude_scalar_dual);
    free(longitude_scalar_dual);
    free(z_scalar_dual);
    free(z_vector_dual);
    free(normal_distance_dual);
    free(f_vec);
    free(to_index);
    free(from_index);
    free(to_index_dual);
    free(from_index_dual);
    free(adjacent_vector_indices_h);
    free(vorticity_indices_pre);
    free(vorticity_indices);
    free(h_curl_indices);
    free(trsk_modified_velocity_indices);
    free(trsk_modified_curl_indices);
    free(adjacent_signs_h);
    free(vorticity_signs_pre);
    free(vorticity_signs);
    free(h_curl_signs);
    free(area_dual_pre);
    free(area_dual);
	free(z_surface);
    return 0;
}






