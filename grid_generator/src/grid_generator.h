/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int find_triangle_on_face_index_from_dual_scalar_on_face_index(int, int, int *, int *, int *, int *);
int find_triangle_edge_points_from_dual_scalar_on_face_index(int, int, int, int *, int *, int *, int [][3], int face_edges[][3], int face_edges_reverse[][3]);
int find_points_per_edge(int, int *);
int find_scalar_points_per_inner_face(int, int *);
int upscale_scalar_point(int, int, int *);
int write_scalar_coordinates(int, int, int, int, int, int, int, double [], double [], double [], double [], double []);
int find_triangles_per_face(int, int *);
int find_v_vector_indices_for_dual_scalar_z(int [], int [], int [], int, int []);
int find_angle_change(double, double, double *);
int find_coords_from_triangle_on_face_index(int, int, int *, int *, int *);
int find_triangle_on_face_index_from_coords(int, int, int, int *);
int find_triangle_indices_from_h_vector_index(int, int, int *, int *, int *, int *, int *, int *, int *, int *, int [][3], int [][3], int [][2], int [][3]);
int find_triangle_edge_points(int, int, int, int *, int *, int *, int *, int *, int *, int *, int [][3], int [][3], int [][3]);
int build_icosahedron(double [], double [], int [][2], int [][3], int [][3], int [][3]);
int generate_horizontal_generators(double [], double [], double [], double [], double [], double [], double [], int [][3], int [][3], int [][3]);
int calc_kinetic_energy(double [], double [], int [], double [], double [], int [], double, int [], int [], double [], double [], double [], double [], int [], double [], double [], double [], double [], int [], int [], double []);
int calc_coriolis_weights(int [], int [], int [], double [], int [], int [], double [], double [], int [], double [], double [], double [], double [], double [], double [], double [], double [], double [], int [], int [], int [], double [], double [], double [], double [], int [], double);
int determine_z_scalar(double [], double [], double [], double, double);
int set_f_vec(double [], double [], double [], double []);











