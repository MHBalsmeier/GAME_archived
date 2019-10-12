enum grid_integers {
NUMBER_OF_PENTAGONS = 12,
RES_ID = 2,
NUMBER_OF_HEXAGONS = (int)(10*(pow(4, RES_ID) - 1)),
NUMBER_OF_EDGES = 30,
NUMBER_OF_SCALARS_H = 12 + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_LAYERS = 2+6*RES_ID,
NUMBER_OF_TRIANGLES = (int)(20*(pow(4, RES_ID))),
NUMBER_OF_SCALARS = NUMBER_OF_SCALARS_H*NUMBER_OF_LAYERS,
NUMBER_OF_VECTORS = NUMBER_OF_SCALARS_H*(NUMBER_OF_LAYERS-1)+NUMBER_OF_VECTORS_H*NUMBER_OF_LAYERS,
NUMBER_OF_DUAL_SCALARS_H = NUMBER_OF_TRIANGLES,
NUMBER_OF_DUAL_SCALARS = NUMBER_OF_DUAL_SCALARS_H*(NUMBER_OF_LAYERS+1),
NUMBER_OF_DUAL_VECTORS_H = 3*(NUMBER_OF_TRIANGLES/2),
NUMBER_OF_DUAL_VECTORS = NUMBER_OF_DUAL_SCALARS_H*NUMBER_OF_LAYERS+NUMBER_OF_DUAL_VECTORS_H*(NUMBER_OF_LAYERS+1),
TRIANGLES_PER_FACE = NUMBER_OF_TRIANGLES/20,
POINTS_PER_EDGE = (int)(pow(2,RES_ID)-1),
POINTS_PER_INNER_FACES = (int)(10*(pow(2,RES_ID)-2)*(pow(2,RES_ID)-1)),
VECTOR_POINTS_PER_INNER_FACE = (int)((pow(2,RES_ID) - 1)*pow(2,RES_ID))};

typedef double Scalar_field[NUMBER_OF_SCALARS];
typedef double Vector_field[NUMBER_OF_VECTORS];
typedef double Dual_scalar_field[NUMBER_OF_DUAL_SCALARS];
typedef double Dual_vector_field[NUMBER_OF_DUAL_VECTORS];

typedef struct grid {
Scalar_field latitude_scalar;
Scalar_field longitude_scalar;
Scalar_field z_scalar;
Vector_field latitude_vector;
Vector_field longitude_vector;
Vector_field z_vector;
Vector_field parallel_distance;
Vector_field normal_distance;
Vector_field gravity;
Scalar_field volume;
Vector_field area;
Vector_field direction[2];
long to_indices[NUMBER_OF_VECTORS];
long from_indices[NUMBER_OF_VECTORS];
long adjacent_vector_indices_h[NUMBER_OF_SCALARS][6];
long adjacent_vector_index_lower[NUMBER_OF_SCALARS];
long adjacent_vector_index_upper[NUMBER_OF_SCALARS];
long adjacent_scalar_indices_h[NUMBER_OF_SCALARS][6];
long adjacent_scalar_index_lower[NUMBER_OF_SCALARS];
long adjacent_scalar_indices_upper[NUMBER_OF_SCALARS];
long adjacent_signs_h[NUMBER_OF_SCALARS][6];
long vorticity_indices[NUMBER_OF_SCALARS][6];
long vorticity_signs[NUMBER_OF_SCALARS][6];
long h_curl_indices[NUMBER_OF_VECTORS][4];
long h_curl_signs[NUMBER_OF_VECTORS][4];
long is_pentagon[NUMBER_OF_SCALARS];
short is_vertical[NUMBER_OF_VECTORS];
long vertical_horizontal_index[NUMBER_OF_VECTORS];
long vector_product_sign[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H];
long recov_hor_par_dual_index[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][11];
double recov_hor_par_dual_weight[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][11];
long recov_hor_ver_dual_index[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][2];
double recov_hor_ver_dual_weight[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][2];
long recov_hor_par_pri_index[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][2];
double recov_hor_par_pri_weight[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][2];
long recov_hor_ver_pri_index[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][4];
double recov_hor_ver_pri_weight[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H][4];
long recov_ver_1_pri_index[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
double recov_ver_1_pri_weight[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
long recov_ver_1_dual_index[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
double recov_ver_1_dual_weight[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
long recov_ver_2_pri_index[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
double recov_ver_2_pri_weight[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
long recov_ver_2_dual_index[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
double recov_ver_2_dual_weight[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H][6];
} Grid;

typedef struct dualgrid {
Dual_scalar_field latitude_scalar;
Dual_scalar_field longitude_scalar;
Dual_scalar_field z_scalar;
Dual_vector_field latitude_vector;
Dual_vector_field longitude_vector;
Dual_vector_field z_vector;
Vector_field direction;
Vector_field normal_distance;
Dual_vector_field area;
Dual_vector_field f_vec;
long to_indices[NUMBER_OF_DUAL_VECTORS];
long from_indices[NUMBER_OF_DUAL_VECTORS];
long adjacent_vector_indices_h[NUMBER_OF_DUAL_SCALARS][3];
long adjacent_vector_index_upper[NUMBER_OF_DUAL_SCALARS];
long adjacent_vector_index_lower[NUMBER_OF_DUAL_SCALARS];
long adjacent_scalar_indices_h[NUMBER_OF_DUAL_SCALARS][3];
long adjacent_scalar_index_lower[NUMBER_OF_DUAL_SCALARS];
long adjacent_scalar_index_upper[NUMBER_OF_DUAL_SCALARS];
long adjacent_signs_h[NUMBER_OF_DUAL_SCALARS][3];
long vorticity_indices[NUMBER_OF_TRIANGLES*NUMBER_OF_LAYERS][3];
long vorticity_signs[NUMBER_OF_TRIANGLES*NUMBER_OF_LAYERS][3];
long h_curl_indices[NUMBER_OF_DUAL_VECTORS_H*(NUMBER_OF_LAYERS+1)][4];
long h_curl_signs[NUMBER_OF_DUAL_VECTORS_H*(NUMBER_OF_LAYERS+1)][4];
short is_vertical[NUMBER_OF_DUAL_VECTORS];
short vertical_horizontal_index[NUMBER_OF_DUAL_VECTORS];
long parallel_distance[NUMBER_OF_DUAL_VECTORS][2];
} Dualgrid;

typedef struct state {
Scalar_field density;
Scalar_field pot_temp;
Scalar_field pressure;
Vector_field wind;
} State;