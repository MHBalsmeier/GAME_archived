/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
The grid generation procedure is manged from this file. Memory allocation and IO is done here, for the rest, functions are called residing in individual files.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include <geos95.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "include.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
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
    int oro_id;
   	oro_id = strtod(argv[1], NULL);
    int n_iterations;
   	n_iterations = strtod(argv[2], NULL);
    int use_scalar_h_file;
   	use_scalar_h_file = strtod(argv[3], NULL);
    int len = strlen(argv[4]);
    char *scalar_h_file = malloc((len + 1)*sizeof(char));
    strcpy(scalar_h_file, argv[4]);
    double stretching_parameter;
   	stretching_parameter = strtof(argv[5], NULL);
   	int no_of_oro_layers = strtod(argv[6], NULL);
   	const int vert_grid_type = strtod(argv[7], NULL);
    double toa;
   	toa = strtof(argv[8], NULL);
   	if (vert_grid_type == 1)
   	{
   		no_of_oro_layers = 0;
   	}
    
    // checking if the no_of_oro_layers is valid
    if (no_of_oro_layers < 0 || no_of_oro_layers >= NO_OF_LAYERS)
    {
    	printf("It must be 0 <= orography_layers < NO_OF_LAYERS.\n");
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
    
	if (no_of_oro_layers >= NO_OF_LAYERS)
	{
		printf("It is no_of_oro_layers >= NO_OF_LAYERS.\n");
		exit(1);
	}
	char grid_name_pre[200];
	char output_file_pre[200];
	char statistics_file_pre[200];
	sprintf(grid_name_pre, "RES%d_L%d_ORO%d", RES_ID, NO_OF_LAYERS, oro_id);
	sprintf(output_file_pre, "grids/RES%d_L%d_ORO%d.nc", RES_ID, NO_OF_LAYERS, oro_id);
	sprintf(statistics_file_pre, "statistics/RES%d_L%d_ORO%d.txt", RES_ID, NO_OF_LAYERS, oro_id);
    char grid_name[strlen(grid_name_pre) + 1];
    char output_file[strlen(output_file_pre) + 1];
    char statistics_file[strlen(statistics_file_pre) + 1];
    strcpy(grid_name, grid_name_pre);
    strcpy(output_file, output_file_pre);
    strcpy(statistics_file, statistics_file_pre);
	printf("Output will be written to file %s.\n", output_file);
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
	double *inner_product_weights = malloc(8*NO_OF_SCALARS*sizeof(double));
    double *density_to_rhombi_weights = malloc(4*NO_OF_VECTORS_H*sizeof(double));
    double *interpol_weights = malloc(5*NO_OF_LATLON_IO_POINTS*sizeof(double));
    double *exner_bg = malloc(NO_OF_SCALARS*sizeof(double));
    double *theta_bg = malloc(NO_OF_SCALARS*sizeof(double));
	double *oro = calloc(NO_OF_SCALARS_H, sizeof(double));
   	double *roughness_length = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sfc_albedo = calloc(NO_OF_SCALARS_H, sizeof(double));
	double *sfc_rho_c = calloc(NO_OF_SCALARS_H, sizeof(double));
	double *t_conductivity = calloc(NO_OF_SCALARS_H, sizeof(double));
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
    int *interpol_indices = malloc(5*NO_OF_LATLON_IO_POINTS*sizeof(int));
	int *is_land = calloc(NO_OF_SCALARS_H, sizeof(int));
    printf(GREEN "finished.\n" RESET);
    
    /*
	1.) creating or reading the properties that determine the horizontal grid
	    ---------------------------------------------------------------------
	*/
    printf("Establishing horizontal grid structure ... \n");
   	int no_of_lloyd_iterations = 0;
    if (use_scalar_h_file == 0)
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
    	read_horizontal_explicit(latitude_scalar, longitude_scalar, from_index, to_index, from_index_dual, to_index_dual, scalar_h_file, &no_of_lloyd_iterations);
    }
    
    /*
    2.) finding the neighbouring vector points of the cells
        ---------------------------------------------------
    */
	find_adjacent_vector_indices_h(from_index, to_index, adjacent_signs_h, adjacent_vector_indices_h);
	
	/*
	3.) grid optimization
	    -----------------
	*/
	if (n_iterations > 0)
	{
		optimize_to_scvt(latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, n_iterations,
		face_edges, face_edges_reverse, face_vertices, adjacent_vector_indices_h, from_index_dual, to_index_dual);
		no_of_lloyd_iterations = no_of_lloyd_iterations + n_iterations;
	}
	
	/*
	4.) determining implicit quantities of the horizontal grid
	    ------------------------------------------------------
	*/
	// calculation of the horizontal coordinates of the dual scalar points
	set_scalar_h_dual_coords(latitude_scalar_dual, longitude_scalar_dual, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices);
	
	// calculation of the horizontal coordinates of the vector points
	set_vector_h_doubles(from_index, to_index, latitude_scalar, longitude_scalar, latitude_vector, longitude_vector, direction);
	
	// Same as before, but for dual vectors. They have the same positions as the primal vectors.
	set_dual_vector_h_doubles(latitude_scalar_dual, latitude_vector, direction_dual, longitude_vector,
	to_index_dual, from_index_dual, longitude_scalar_dual, rel_on_line_dual);
	
	// determining the directions of the dual vectors
	direct_tangential_unity(latitude_scalar_dual, longitude_scalar_dual, direction, direction_dual,
	to_index_dual, from_index_dual, rel_on_line_dual, ORTH_CRITERION_DEG);
	
	// setting the Coriolis vector
    set_f_vec(latitude_vector, direction, direction_dual, f_vec);
    
    // calculating the dual cells on the unity sphere
    calc_triangle_area_unity(triangle_face_unit_sphere, latitude_scalar, longitude_scalar, face_edges,
    face_edges_reverse, face_vertices);
    
    // finding the vorticity indices
	calc_vorticity_indices_triangles(from_index_dual, to_index_dual, direction, direction_dual,
	vorticity_indices_triangles, ORTH_CRITERION_DEG, vorticity_signs_triangles);
	
	// calculating the cell faces on the unity sphere
	calc_cell_area_unity(pent_hex_face_unity_sphere, latitude_scalar_dual,
	longitude_scalar_dual, adjacent_vector_indices_h, vorticity_indices_triangles);
    printf(GREEN "Horizontal grid structure determined.\n" RESET);
	
	/*
	5.) setting the physical surface properties
	    ---------------------------------------
	*/
	set_sfc_properties(latitude_scalar, longitude_scalar, roughness_length, sfc_albedo, sfc_rho_c, t_conductivity, oro, is_land, oro_id);
	
	/*
	6.) setting the explicit property of the vertical grid
	    --------------------------------------------------
	*/
    printf("Setting the vertical coordinates of the scalar data points ... ");
	set_z_scalar(z_scalar, oro, no_of_oro_layers, toa, stretching_parameter, vert_grid_type);
    printf(GREEN "finished.\n" RESET);
	
	/*
	7.) setting the implicit quantities of the vertical grid
	    ----------------------------------------------------
	*/
	if (vert_grid_type == 1)
	{
		set_scalar_shading_indices(z_scalar, oro, no_of_shaded_points_scalar);
		set_vector_shading_indices(from_index, to_index, no_of_shaded_points_scalar, no_of_shaded_points_vector);
	}
	
	printf("Determining vector z coordinates and normal distances of the primal grid ... ");
	set_z_vector_and_normal_distance(z_vector, z_scalar, normal_distance, latitude_scalar, longitude_scalar, from_index, to_index, toa, vert_grid_type, oro);
	free(oro);
    printf(GREEN "finished.\n" RESET);	
	
	printf("Determining scalar z coordinates of the dual grid ... ");
	set_z_scalar_dual(z_scalar_dual, z_vector, from_index, to_index, vorticity_indices_triangles, toa);
    printf(GREEN "finished.\n" RESET);
	
	printf("Determining vector z coordinates of the dual grid and distances of the dual grid ... ");
	calc_z_vector_dual_and_normal_distance_dual(z_vector_dual, normal_distance_dual, z_scalar_dual, toa, from_index, to_index, z_vector,
	from_index_dual, to_index_dual, latitude_scalar_dual, longitude_scalar_dual, vorticity_indices_triangles);
    printf(GREEN "finished.\n" RESET);
	
	printf("Calculating areas ... ");
	set_area(area, z_vector, z_vector_dual, normal_distance_dual, pent_hex_face_unity_sphere);
    printf(GREEN "finished.\n" RESET);
    
    printf("Calculating dual areas ... ");
	set_area_dual(area_dual, z_vector_dual, normal_distance, z_vector, from_index, to_index, triangle_face_unit_sphere, toa);
    printf(GREEN "finished.\n" RESET);
    
    printf("Calculating grid box volumes ... ");
	set_volume(volume, z_vector, area, from_index, to_index, toa, vorticity_indices_triangles);
    printf(GREEN "finished.\n" RESET);
    
    /*
    8.) Now come the derived quantities, which are needed for differential operators.
        -----------------------------------------------------------------------------
    */
    printf("Setting the gravity potential ... ");
	set_gravity_potential(z_scalar, gravity_potential, GRAVITY_MEAN_SFC_ABS);
    printf(GREEN "finished.\n" RESET);
    
    printf("Setting the hydrostatic background state ... ");
	set_background_state(z_scalar, gravity_potential, theta_bg, exner_bg);
    printf(GREEN "finished.\n" RESET);
    
    printf("Calculating inner product weights ... ");
	calc_inner_product(inner_product_weights, normal_distance, volume, to_index, from_index, area, z_scalar, z_vector, adjacent_vector_indices_h);
    printf(GREEN "finished.\n" RESET);
    
    printf("Setting rhombus interpolation indices and weights ... ");
	rhombus_averaging(vorticity_indices_triangles, vorticity_signs_triangles, from_index_dual,
	to_index_dual, vorticity_indices_rhombi, density_to_rhombi_indices, from_index, to_index, area_dual,
	z_vector, latitude_scalar_dual, longitude_scalar_dual, density_to_rhombi_weights, latitude_vector,
	longitude_vector, latitude_scalar, longitude_scalar, toa);
    printf(GREEN "finished.\n" RESET);
    
    printf("Calculating Coriolis indices and weights ... ");
	coriolis(from_index_dual, to_index_dual, trsk_modified_curl_indices, normal_distance, normal_distance_dual,
	to_index, area, z_scalar, latitude_scalar, longitude_scalar, latitude_vector, longitude_vector, latitude_scalar_dual,
	longitude_scalar_dual, trsk_weights, trsk_indices, from_index, adjacent_vector_indices_h, z_vector, z_vector_dual);
    printf(GREEN "finished.\n" RESET);
    
    printf("Calculating interpolation to the lat-lon grid ... ");
    interpolate_ll(latitude_scalar, longitude_scalar, interpol_indices, interpol_weights);
    printf(GREEN "finished.\n" RESET);
    
    // A statistics file is created to compare the fundamental statistical properties of the grid with the literature.
	write_statistics_file(pent_hex_face_unity_sphere, normal_distance, normal_distance_dual, no_of_lloyd_iterations, grid_name, statistics_file);
	
	/*
	writing the result to a netcdf file
    */
    int retval, latitude_scalar_id, longitude_scalar_id, direction_id, latitude_vector_id, longitude_vector_id, latitude_scalar_dual_id, longitude_scalar_dual_id,
    z_scalar_id, z_vector_id, normal_distance_id, volume_id, area_id, trsk_weights_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id,
    from_index_id, to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, trsk_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id,
    vorticity_signs_triangles_id, f_vec_dimid, scalar_dimid, scalar_h_dimid, scalar_dual_h_dimid, vector_dimid, latlon_dimid_5, scalar_h_dimid_6, vector_h_dimid,
    vector_h_dimid_10, vector_h_dimid_4, vector_v_dimid_6, vector_dual_dimid, gravity_potential_id, scalar_dual_h_dimid_3, vector_dual_area_dimid,
    inner_product_weights_id, scalar_8_dimid, scalar_2_dimid, vector_h_dual_dimid_2, density_to_rhombi_indices_id, density_to_rhombi_weights_id,
    vorticity_indices_triangles_id, ncid_g_prop, single_double_dimid, no_of_shaded_points_vector_id, no_of_shaded_points_scalar_id,
    no_of_lloyd_iterations_id, single_int_dimid, interpol_indices_id, interpol_weights_id, theta_bg_id, exner_bg_id, sfc_albedo_id, sfc_rho_c_id, t_conductivity_id,
    roughness_length_id, is_land_id, vert_grid_type_id, no_of_oro_layers_id, stretching_parameter_id, toa_id;
    
    printf("Starting to write to output file ... ");
    if ((retval = nc_create(output_file, NC_CLOBBER, &ncid_g_prop)))
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
    if ((retval = nc_def_dim(ncid_g_prop, "latlon_3_index", 5*NO_OF_LATLON_IO_POINTS, &latlon_dimid_5)))
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
    if ((retval = nc_def_var(ncid_g_prop, "no_of_lloyd_iterations", NC_INT, 1, &single_int_dimid, &no_of_lloyd_iterations_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vert_grid_type", NC_INT, 1, &single_int_dimid, &vert_grid_type_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "no_of_oro_layers", NC_INT, 1, &single_int_dimid, &no_of_oro_layers_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "stretching_parameter", NC_DOUBLE, 1, &single_double_dimid, &stretching_parameter_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "toa", NC_DOUBLE, 1, &single_double_dimid, &toa_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "theta_bg", NC_DOUBLE, 1, &scalar_dimid, &theta_bg_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, theta_bg_id, "units", strlen("K"), "K")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "exner_bg", NC_DOUBLE, 1, &scalar_dimid, &exner_bg_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "interpol_weights", NC_DOUBLE, 1, &latlon_dimid_5, &interpol_weights_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "interpol_indices", NC_INT, 1, &latlon_dimid_5, &interpol_indices_id)))
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
	if ((retval = nc_def_var(ncid_g_prop, "sfc_albedo", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_albedo_id)))
	  	ERR(retval);
	if ((retval = nc_def_var(ncid_g_prop, "sfc_rho_c", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_rho_c_id)))
	  	ERR(retval);
	if ((retval = nc_put_att_text(ncid_g_prop, sfc_rho_c_id, "units", strlen("J/(K*m**3)"), "J/(K*m**3)")))
	  	ERR(retval);
	if ((retval = nc_def_var(ncid_g_prop, "is_land", NC_INT, 1, &scalar_h_dimid, &is_land_id)))
	  	ERR(retval);
	if ((retval = nc_def_var(ncid_g_prop, "t_conductivity", NC_DOUBLE, 1, &scalar_h_dimid, &t_conductivity_id)))
	  	ERR(retval);
	if ((retval = nc_put_att_text(ncid_g_prop, t_conductivity_id, "units", strlen("m^2/2"), "m^2/2")))
	  	ERR(retval);
	if ((retval = nc_def_var(ncid_g_prop, "roughness_length", NC_DOUBLE, 1, &scalar_h_dimid, &roughness_length_id)))
	  	ERR(retval);
	if ((retval = nc_put_att_text(ncid_g_prop, roughness_length_id, "units", strlen("m"), "m")))
	  	ERR(retval);
    if ((retval = nc_enddef(ncid_g_prop)))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_shaded_points_scalar_id, &no_of_shaded_points_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_shaded_points_vector_id, &no_of_shaded_points_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_lloyd_iterations_id, &no_of_lloyd_iterations)))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vert_grid_type_id, &vert_grid_type)))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, stretching_parameter_id, &stretching_parameter)))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, toa_id, &toa)))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, no_of_oro_layers_id, &no_of_oro_layers)))
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
    if ((retval = nc_put_var_double(ncid_g_prop, theta_bg_id, &theta_bg[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, exner_bg_id, &exner_bg[0])))
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
	if ((retval = nc_put_var_double(ncid_g_prop, sfc_albedo_id, &sfc_albedo[0])))
	  	ERR(retval);
	if ((retval = nc_put_var_double(ncid_g_prop, sfc_rho_c_id, &sfc_rho_c[0])))
	  	ERR(retval);
	if ((retval = nc_put_var_double(ncid_g_prop, t_conductivity_id, &t_conductivity[0])))
	  	ERR(retval);
	if ((retval = nc_put_var_double(ncid_g_prop, roughness_length_id, &roughness_length[0])))
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
	if ((retval = nc_put_var_int(ncid_g_prop, is_land_id, &is_land[0])))
	  	ERR(retval);
    if ((retval = nc_close(ncid_g_prop)))
        ERR(retval);
    printf(GREEN "finished.\n" RESET);
    
    // freeing allocated memory
    free(scalar_h_file);
	free(roughness_length);
	free(sfc_albedo);
	free(sfc_rho_c);
	free(t_conductivity);
	free(is_land);
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
    free(exner_bg);
    free(theta_bg);
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
	free(interpol_indices);
	free(interpol_weights);
    return 0;
}






