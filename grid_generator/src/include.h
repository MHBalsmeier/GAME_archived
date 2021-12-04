/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int find_triangle_on_face_index_from_dual_scalar_on_face_index(int, int, int *, int *, int *, int *);
int find_triangle_edge_points_from_dual_scalar_on_face_index(int, int, int, int *, int *, int *, int [][3], int face_edges[][3], int face_edges_reverse[][3]);
int find_points_per_edge(int, int *);
int find_scalar_points_per_inner_face(int, int *);
int upscale_scalar_point(int, int, int *);
int write_scalar_coordinates(int, int, int, int, int, int, int, double [], double [], double [], double [], double []);
int find_triangles_per_face(int, int *);
int find_v_vector_indices_for_dual_scalar_z(int [], int [], int [], int, int []);
int find_coords_from_triangle_on_face_index(int, int, int *, int *, int *);
int find_triangle_on_face_index_from_coords(int, int, int, int *);
int find_triangle_indices_from_h_vector_index(int, int, int *, int *, int *, int *, int *, int *, int *, int *, int [][3], int [][3], int [][2], int [][3]);
int find_triangle_edge_points(int, int, int, int *, int *, int *, int *, int *, int *, int *, int [][3], int [][3], int [][3]);
int build_icosahedron(double [], double [], int [][2], int [][3], int [][3], int [][3]);
int generate_horizontal_generators(double [], double [], double [], double [], double [], double [], double [], int [][3], int [][3], int [][3]);
int calc_inner_product(double [], double [], double [], int [], int [], double [], double [], double [], int []);
int coriolis(int [], int [], int [], double [], double [], int [], double [], double [], double [], double [], double [], double [], double [], double [], double [], int [], int [], int [], double [], double []);
int set_z_scalar(double [], double [], int, double, double, int);
int set_f_vec(double [], double [], double [], double []);
int set_orography(int, int, double []);
int calc_cell_area_unity(double [], double [], double [], int [], int []);
int calc_triangle_area_unity(double [], double [], double [], int [][3], int [][3], int [][3], int [][2]);
int set_vector_h_doubles(int [], int [], double [], double [], double [], double [], double []);
int set_from_to_index(int [], int [], int [][3], int [][3], int [][3], int [][2]);
int set_scalar_h_dual_coords(double [], double [], double [], double [], int [][3], int [][3], int [][3], int [][2]);
int set_z_scalar_dual(double [], double [], int [], int [], int [], double);
int set_volume(double [], double [], double [], int [], int [], double, int []);
int find_adjacent_vector_indices_h(int [], int [], int [], int []);
int set_horizontal_curl_indices(double [], double [], int [], int [], int [], double, int []);
int rhombus_averaging(int [], int [], int [], int [], int [], int [], int [], int [], double [], double [], double [], double [], double [], double [], double [], double [], double [], double);
int set_dual_vector_h_doubles(double [], double [], double [], double [], int [], int [], double [], double []);
int set_gravity_potential(double [], double [], double);
int set_from_to_index_dual(int [], int [], int [][3], int [][3]);
int calc_vorticity_indices_triangles(int [], int [], double [], double [], int [], double, int []);
int set_z_vector_and_normal_distance(double [], double [], double [], double [], double [], int [], int [], double, int, double []);
int set_area(double [], double [], double [], double [], double []);
int calc_z_vector_dual_and_normal_distance_dual(double [], double [], double [], double, int [], int [], double [], int [], int [], double [], double [], int []);
int set_area_dual(double [], double [], double [], double [], int [], int [], double [], double);
int optimize_to_scvt(double [], double [], double [], double [], int, int [][3], int [][3], int [][3], int [][2], int [], int [], int []);
int read_horizontal_explicit(double [], double [], int [], int [], int [], int [], char [], int *);
int write_statistics_file(double [], double [], double [], char []);
int set_scalar_shading_indices(double [], double [], int []);
int set_vector_shading_indices(int [], int [], int [], int []);
int direct_tangential_unity(double [], double [], double [], double [], int [], int [], double [], double);
int interpolate_ll(double [], double [], int [], double []);
int set_background_state(double [], double [], double [], double []);











