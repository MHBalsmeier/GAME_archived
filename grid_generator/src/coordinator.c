/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
The grid generation procedure is manged from this file. Memory allocation and IO is done here, for the rest, functions are called residing in individual files.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include "include.h"
#include "enum.h"
#include "geos95.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
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
    double stretching_parameter;
   	stretching_parameter = strtof(argv[6], NULL);
   	int NO_OF_ORO_LAYERS = strtod(argv[7], NULL);
   	const int VERT_GRID_TYPE = strtod(argv[8], NULL);
    double TOA;
   	TOA = strtof(argv[9], NULL);
   	if (VERT_GRID_TYPE == 1)
   	{
   		NO_OF_ORO_LAYERS = 0;
   	}
    
    // Checking wether the RES_ID of the SCALAR_H_FILE corresponds to the RES_ID in enum.h.
    char res_id_as_string[2];
    res_id_as_string[0] = SCALAR_H_FILE[7];
    res_id_as_string[1] = SCALAR_H_FILE[8];
    int res_id_from_scalar_h_file;
    if (res_id_as_string[1] == 'L')
    {
    	res_id_from_scalar_h_file = res_id_as_string[0] - '0';
    }
    else
    {
    	res_id_from_scalar_h_file = strtod(res_id_as_string, NULL);
    }
    if (USE_SCALAR_H_FILE == 1 && res_id_from_scalar_h_file != RES_ID)
    {
    	printf("The resolution (res_id = %d) of the scalar_h_coords_file does not correspond to the resolution (RES_ID = %d) in enum.h.\n", res_id_from_scalar_h_file, RES_ID);
    	printf("Recompile with RES_ID = %d or choose a scalar_h_coords_file with res_id = %d, then try again.\n", res_id_from_scalar_h_file, RES_ID);
    	printf("Aborting.\n");
    	exit(1);
    }
    
    // cechking wether the stretching parameter is in a valid range
    if (stretching_parameter < 1)
    {
    	printf("stretching_parameter must be >= 1.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    
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
    double *trsk_weights = calloc(10*NO_OF_VECTORS_H, sizeof(double));
    double *latitude_scalar_dual = malloc(NO_OF_DUAL_SCALARS_H*sizeof(double));
    double *longitude_scalar_dual = malloc(NO_OF_DUAL_SCALARS_H*sizeof(double));
    double *z_scalar_dual = malloc(NO_OF_DUAL_SCALARS*sizeof(double));
    double *z_vector_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *normal_distance_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *direction_dual = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *area_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *f_vec = malloc(2*NO_OF_VECTORS_H*sizeof(double));
    double *triangle_face_unit_sphere = malloc(NO_OF_DUAL_SCALARS_H*sizeof(double));
    double *pent_hex_face_unity_sphere = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *rel_on_line_dual = malloc(NO_OF_VECTORS_H*sizeof(double));
	double *z_surface = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *inner_product_weights = malloc(8*NO_OF_SCALARS*sizeof(double));
    double *density_to_rhombi_weights = malloc(4*NO_OF_VECTORS_H*sizeof(double));
    double *interpol_weights = malloc(3*NO_OF_LATLON_IO_POINTS*sizeof(double));
    int *to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *from_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *trsk_indices = calloc(10*NO_OF_VECTORS_H, sizeof(int));
    int *trsk_modified_curl_indices = calloc(10*NO_OF_VECTORS_H, sizeof(int));
    int *no_of_shaded_points_scalar = calloc(NO_OF_SCALARS_H, sizeof(int));
    int *no_of_shaded_points_vector = calloc(NO_OF_VECTORS_H, sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *vorticity_indices_triangles = malloc(3*NO_OF_DUAL_SCALARS_H*sizeof(int));
    int *vorticity_indices_rhombi = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *to_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *from_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_signs_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *vorticity_signs_triangles = malloc(3*NO_OF_DUAL_SCALARS_H*sizeof(int));
    int *density_to_rhombi_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *interpol_indices = malloc(3*NO_OF_LATLON_IO_POINTS*sizeof(int));
    printf(GREEN "finished.\n" RESET);
    
    /*
    1.) reading the orography
        ---------------------
    */
    printf("Reading orography data ... ");
	set_orography(RES_ID, ORO_ID, z_surface);
    printf(GREEN "finished.\n" RESET);
    
    /*
	2.) creating or reading the properties that determine the horizontal grid
	    ---------------------------------------------------------------------
	*/
    printf("Establishing horizontal grid structure ... \n");
   	int no_of_lloyd_cycles = 0;
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
    	read_horizontal_explicit(latitude_scalar, longitude_scalar, from_index, to_index, from_index_dual, to_index_dual, SCALAR_H_FILE, &no_of_lloyd_cycles);
    }
    
    /*
    3.) finding the neighbouring vector points of the cells
        ---------------------------------------------------
    */
	find_adjacent_vector_indices_h(from_index, to_index, adjacent_signs_h, adjacent_vector_indices_h);
	
	/*
	4.) grid optimization
	    -----------------
	*/
	if (OPTIMIZE_BOOL == 1)
	{
		optimize_to_scvt(latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, N_ITERATIONS, face_edges, face_edges_reverse, face_vertices, edge_vertices, adjacent_vector_indices_h, from_index_dual, to_index_dual);
		no_of_lloyd_cycles = no_of_lloyd_cycles + N_ITERATIONS;
	}
	
	/*
	5.) determining implicit quantities of the horizontal grid
	    ------------------------------------------------------
	*/
	// Calculation of the horizontal coordinates of the dual scalar points. The dual scalar points are the vertices of the Voronoi mesh of the primal grid.
	set_scalar_h_dual_coords(latitude_scalar_dual, longitude_scalar_dual, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices, edge_vertices);
	// Calculation of the horizontal coordinates of the vector points. The vector points are in the middle between the neighbouring scalar points.
	set_vector_h_doubles(from_index, to_index, latitude_scalar, longitude_scalar, latitude_vector, longitude_vector, direction);
	// Same as before but for dual vectors. They have the same positions as the primal vectors.
	set_dual_vector_h_doubles(latitude_scalar_dual, latitude_vector, direction_dual, longitude_vector,
	to_index_dual, from_index_dual, longitude_scalar_dual, rel_on_line_dual);
	direct_tangential_unity(latitude_scalar_dual, longitude_scalar_dual, direction, direction_dual,
	to_index_dual, from_index_dual, rel_on_line_dual, ORTH_CRITERION_DEG);
	// Setting the Coriolis vector.
    set_f_vec(latitude_vector, direction, direction_dual, f_vec);
    // calculating the dual cells on the unity sphere
    calc_triangle_face_unity(triangle_face_unit_sphere, latitude_scalar, longitude_scalar, face_edges,
    face_edges_reverse, face_vertices, edge_vertices);
    // finding the vorticity indices
	calc_vorticity_indices_triangles(from_index_dual, to_index_dual, direction, direction_dual,
	vorticity_indices_triangles, ORTH_CRITERION_DEG, vorticity_signs_triangles);
	// checking if the grid is orthogonal
	check_for_orthogonality(direction, direction_dual, ORTH_CRITERION_DEG);
	// calculating the cell faces on the unity sphere
	calc_cell_face_unity(pent_hex_face_unity_sphere, latitude_scalar_dual,
	longitude_scalar_dual, adjacent_vector_indices_h, vorticity_indices_triangles);
    printf(GREEN "Horizontal grid structure determined.\n" RESET);
	
	/*
	6.) setting the explicit property of the vertical grid
	    --------------------------------------------------
	*/
    printf("Setting the vertical coordinates of the scalar data points ... ");
	set_z_scalar(z_scalar, z_surface, NO_OF_ORO_LAYERS, TOA, stretching_parameter, VERT_GRID_TYPE);
    printf(GREEN "finished.\n" RESET);
	
	/*
	7.) setting the implicit quantities of the vertical grid
	    ----------------------------------------------------
	*/
    // setting the gravity potential
    printf("Setting gravity potential ... ");
	set_gravity_potential(z_scalar, gravity_potential, GRAVITY_MEAN_SFC_ABS);
    printf(GREEN "finished.\n" RESET);
	if (VERT_GRID_TYPE == 1)
	{
		set_scalar_shading_indices(z_scalar, z_surface, no_of_shaded_points_scalar);
		set_vector_shading_indices(from_index, to_index, no_of_shaded_points_scalar, no_of_shaded_points_vector);
	}
	set_z_vector_and_normal_distance(z_vector, z_scalar, normal_distance, latitude_scalar, longitude_scalar, from_index, to_index, TOA, VERT_GRID_TYPE, z_surface);
	printf("Mapping horizontal areas from unit sphere to model levels ... ");
	map_hor_area_to_half_levels(area, z_vector, pent_hex_face_unity_sphere);
    printf(GREEN "finished.\n" RESET);
	printf("Determining scalar z coordinates of the dual grid ... ");
	set_z_scalar_dual(z_scalar_dual, z_vector, from_index, to_index, vorticity_indices_triangles, TOA);
    printf(GREEN "finished.\n" RESET);
    printf("Calculating grid box volumes ... ");
	set_volume(volume, z_vector, area, from_index, to_index, TOA, vorticity_indices_triangles);
    printf(GREEN "finished.\n" RESET);
	printf("Determining vector z coordinates of the dual grid and distances of the dual grid ... ");
	calc_z_vector_dual_and_normal_distance_dual(z_vector_dual, normal_distance_dual, z_scalar_dual, TOA, from_index, to_index, z_vector,
	from_index_dual, to_index_dual, latitude_scalar_dual, longitude_scalar_dual, vorticity_indices_triangles);
    printf(GREEN "finished.\n" RESET);
    printf("Calculating dual areas ... ");
	set_area_dual(area_dual, z_vector_dual, normal_distance, z_vector, from_index, to_index, triangle_face_unit_sphere, TOA);
    printf(GREEN "finished.\n" RESET);
    printf("Calculating vertical faces, pre version ... ");
	calculate_vertical_faces(area, z_vector_dual, normal_distance_dual, TOA);
    printf(GREEN "finished.\n" RESET);
    
    /*
    8.) Now come the derived quantities, which are needed for differential operators.
        -----------------------------------------------------------------------------
    */
    // interpolation to the lat-lon grid
    printf("Calculating interpolation to the lat-lon grid ... ");
    interpolate_ll(latitude_scalar, longitude_scalar, interpol_indices, interpol_weights);
    printf(GREEN "finished.\n" RESET);
    // inner product
    printf("Calculating inner product weights ... ");
	calc_inner_product(inner_product_weights, normal_distance, volume, to_index, from_index, area, z_scalar, z_vector, adjacent_vector_indices_h);
    printf(GREEN "finished.\n" RESET);
    // setting indices and weights related to averaging of a scalar quantity to rhombi
    printf("Setting rhombus interpolation indices and weights ... ");
	rhombus_averaging(vorticity_indices_triangles, vorticity_signs_triangles, from_index_dual,
	to_index_dual, vorticity_indices_rhombi, density_to_rhombi_indices, from_index, to_index, area_dual,
	z_vector, latitude_scalar_dual, longitude_scalar_dual, density_to_rhombi_weights, latitude_vector,
	longitude_vector, latitude_scalar, longitude_scalar, TOA);
    printf(GREEN "finished.\n" RESET);
    // modified TRSK
    printf("Calculating Coriolis indices and weights ... ");
	coriolis(from_index_dual, to_index_dual, trsk_modified_curl_indices, normal_distance, normal_distance_dual,
	to_index, area, z_scalar, latitude_scalar, longitude_scalar, latitude_vector, longitude_vector, latitude_scalar_dual,
	longitude_scalar_dual, trsk_weights, trsk_indices, from_index, adjacent_vector_indices_h, z_vector, z_vector_dual);
    printf(GREEN "finished.\n" RESET);
    
    // A statistics file is created to compare the fundamental statistical properties of the grid with the literature.
	write_statistics_file(pent_hex_face_unity_sphere, normal_distance, normal_distance_dual, STATISTICS_FILE);
	
	/*
	writing the result to a netcdf file
    */
    int retval, latitude_scalar_id, longitude_scalar_id, direction_id, latitude_vector_id, longitude_vector_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_id, z_vector_id, normal_distance_id, volume_id, area_id, trsk_weights_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, trsk_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id, vorticity_signs_triangles_id, f_vec_dimid, scalar_dimid, scalar_h_dimid, scalar_dual_h_dimid, vector_dimid, latlon_dimid_3, scalar_h_dimid_6, vector_h_dimid, vector_h_dimid_10, vector_h_dimid_4, vector_v_dimid_6, vector_dual_dimid, gravity_potential_id, scalar_dual_h_dimid_3, vector_dual_area_dimid, inner_product_weights_id, scalar_8_dimid, scalar_2_dimid, vector_h_dual_dimid_2, density_to_rhombi_indices_id, density_to_rhombi_weights_id, vorticity_indices_triangles_id, ncid_g_prop, single_double_dimid, stretching_parameter_id, no_of_shaded_points_vector_id, no_of_shaded_points_scalar_id, no_of_lloyd_cycles_id, single_int_dimid, interpol_indices_id, interpol_weights_id;
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
    if ((retval = nc_def_dim(ncid_g_prop, "latlon_3_index", 3*NO_OF_LATLON_IO_POINTS, &latlon_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_6_index", 6*NO_OF_SCALARS_H, &scalar_h_dimid_6)))
        ERR(retval);
	if ((retval = nc_def_dim(ncid_g_prop, "vector_h_10_index", 10*NO_OF_VECTORS_H, &vector_h_dimid_10)))
	    ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_4_index", 4*NO_OF_VECTORS_H, &vector_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_6_index", 6*NO_OF_LEVELS*NO_OF_SCALARS_H, &vector_v_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "f_vec_index", 2*NO_OF_VECTORS_H, &f_vec_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual", NO_OF_DUAL_VECTORS, &vector_dual_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual_area", NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS, &vector_dual_area_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_2_dual", 2*NO_OF_DUAL_H_VECTORS, &vector_h_dual_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "single_double_dimid_index", 1, &single_double_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "single_int_dimid_index", 1, &single_int_dimid)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "stretching_parameter", NC_DOUBLE, 1, &single_double_dimid, &stretching_parameter_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "no_of_lloyd_cycles", NC_INT, 1, &single_int_dimid, &no_of_lloyd_cycles_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "z_vector", NC_DOUBLE, 1, &vector_dimid, &z_vector_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "trsk_weights", NC_DOUBLE, 1, &vector_h_dimid_10, &trsk_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector_dual", NC_DOUBLE, 1, &vector_dual_dimid, &z_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance_dual", NC_DOUBLE, 1, &vector_dual_dimid, &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "area_dual", NC_DOUBLE, 1, &vector_dual_dimid, &area_dual_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "inner_product_weights", NC_DOUBLE, 1, &scalar_8_dimid, &inner_product_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "density_to_rhombi_weights", NC_DOUBLE, 1, &vector_h_dimid_4, &density_to_rhombi_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "interpol_weights", NC_DOUBLE, 1, &latlon_dimid_3, &interpol_weights_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "interpol_indices", NC_INT, 1, &latlon_dimid_3, &interpol_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_indices", NC_INT, 1, &vector_h_dimid_10, &trsk_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_curl_indices", NC_INT, 1, &vector_h_dimid_10, &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_signs_h", NC_INT, 1, &scalar_h_dimid_6, &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs_triangles", NC_INT, 1, &scalar_dual_h_dimid_3, &vorticity_signs_triangles_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices_triangles", NC_INT, 1, &scalar_dual_h_dimid_3, &vorticity_indices_triangles_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "density_to_rhombi_indices", NC_INT, 1, &vector_h_dimid_4, &density_to_rhombi_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "no_of_shaded_points_scalar", NC_INT, 1, &scalar_h_dimid, &no_of_shaded_points_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "no_of_shaded_points_vector", NC_INT, 1, &vector_h_dimid, &no_of_shaded_points_vector_id)))
        ERR(retval);
    if ((retval = nc_enddef(ncid_g_prop)))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_shaded_points_scalar_id, &no_of_shaded_points_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_shaded_points_vector_id, &no_of_shaded_points_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_lloyd_cycles_id, &no_of_lloyd_cycles)))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, stretching_parameter_id, &stretching_parameter)))
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
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, volume_id, &volume[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_id, &area[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, inner_product_weights_id, &inner_product_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, trsk_weights_id, &trsk_weights[0])))
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
    if ((retval = nc_put_var_double(ncid_g_prop, density_to_rhombi_weights_id, &density_to_rhombi_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, interpol_weights_id, &interpol_weights[0])))
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
    if ((retval = nc_put_var_int(ncid_g_prop, trsk_indices_id, &trsk_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, trsk_modified_curl_indices_id, &trsk_modified_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vorticity_signs_triangles_id, &vorticity_signs_triangles[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vorticity_indices_triangles_id, &vorticity_indices_triangles[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, density_to_rhombi_indices_id, &density_to_rhombi_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, interpol_indices_id, &interpol_indices[0])))
        ERR(retval);
    if ((retval = nc_close(ncid_g_prop)))
        ERR(retval);
    printf(GREEN "finished.\n" RESET);
    
    // freeing allocated memory
    free(SCALAR_H_FILE);
    free(no_of_shaded_points_scalar);
    free(no_of_shaded_points_vector);
    free(latitude_ico);
    free(longitude_ico);
    free(x_unity);
    free(y_unity);
    free(z_unity);
    free(pent_hex_face_unity_sphere);
    free(triangle_face_unit_sphere);
    free(direction_dual);
    free(density_to_rhombi_weights);
    free(density_to_rhombi_indices);
    free(rel_on_line_dual);
    free(inner_product_weights);
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
    free(trsk_weights);
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
    free(vorticity_indices_triangles);
    free(vorticity_indices_rhombi);
    free(trsk_indices);
    free(trsk_modified_curl_indices);
    free(adjacent_signs_h);
    free(vorticity_signs_triangles);
    free(area_dual);
	free(z_surface);
	free(interpol_indices);
	free(interpol_weights);
    return 0;
}






