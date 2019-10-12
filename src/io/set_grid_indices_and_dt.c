#include <math.h>
#include "../enum_and_typedefs.h"
#include "io.h"

void find_min_dist(void);
double calculate_distance_h(double,double,double,double,double);
double calculate_vertical_face(double, double, double);
void find_midpoint(double,double,double,double,double*,double*);
void find_geos(double,double,double,double*,double*);
void find_geodetic(double,double,double,double,double,double*,double*);
void find_center(double,double,double,double,double,double,double*,double*);
double find_geodetic_direction(double,double,double,double,double);
extern double semimajor;
extern double omega;
extern double min_dist;
extern Grid grid;
extern Dualgrid dualgrid;

void set_grid_properties(void)
{
	double latitude_ico[12];
	latitude_ico[0] = M_PI/2;
	latitude_ico[1] = M_PI/6;
	latitude_ico[2] = M_PI/6;
	latitude_ico[3] = M_PI/6;
	latitude_ico[4] = M_PI/6;
	latitude_ico[5] = M_PI/6;
	latitude_ico[6] = -M_PI/6;
	latitude_ico[7] = -M_PI/6;
	latitude_ico[8] = -M_PI/6;
	latitude_ico[9] = -M_PI/6;
	latitude_ico[10] = -M_PI/6;
	latitude_ico[11] = -M_PI/2;
	double longitude_ico[12];
	longitude_ico[0] = 0;
	longitude_ico[1] = 0;
	longitude_ico[2] = 1*2*M_PI/5;
	longitude_ico[3] = 2*2*M_PI/5;
	longitude_ico[4] = 3*2*M_PI/5;
	longitude_ico[5] = 4*2*M_PI/5;
	longitude_ico[6] = 0;
	longitude_ico[7] = 1*2*M_PI/5+2*M_PI/10;
	longitude_ico[8] = 2*2*M_PI/5+2*M_PI/10;
	longitude_ico[9] = 3*2*M_PI/5+2*M_PI/10;
	longitude_ico[10] = 4*2*M_PI/5+2*M_PI/10;
	longitude_ico[11] = 0;
	double edges_vertices[NUMBER_OF_EDGES][2];
	edges_vertices[0][0] = 0;
	edges_vertices[0][1] = 1;
	edges_vertices[1][0] = 0;
	edges_vertices[1][1] = 2;
	edges_vertices[2][0] = 0;
	edges_vertices[2][1] = 3;
	edges_vertices[3][0] = 0;
	edges_vertices[3][1] = 4;
	edges_vertices[4][0] = 0;
	edges_vertices[4][1] = 5;
	edges_vertices[5][0] = 1;
	edges_vertices[5][1] = 2;
	edges_vertices[6][0] = 2;
	edges_vertices[6][1] = 3;
	edges_vertices[7][0] = 3;
	edges_vertices[7][1] = 4;
	edges_vertices[8][0] = 4;
	edges_vertices[8][1] = 5;
	edges_vertices[9][0] = 5;
	edges_vertices[9][1] = 1;
	edges_vertices[10][0] = 1;
	edges_vertices[10][1] = 6;
	edges_vertices[11][0] = 6;
	edges_vertices[11][1] = 2;
	edges_vertices[12][0] = 2;
	edges_vertices[12][1] = 7;
	edges_vertices[13][0] = 7;
	edges_vertices[13][1] = 3;
	edges_vertices[14][0] = 3;
	edges_vertices[14][1] = 8;
	edges_vertices[15][0] = 8;
	edges_vertices[15][1] = 4;
	edges_vertices[16][0] = 4;
	edges_vertices[16][1] = 9;
	edges_vertices[17][0] = 9;
	edges_vertices[17][1] = 5;
	edges_vertices[18][0] = 5;
	edges_vertices[18][1] = 10;
	edges_vertices[19][0] = 10;
	edges_vertices[19][1] = 1;
	edges_vertices[20][0] = 10;
	edges_vertices[20][1] = 6;
	edges_vertices[21][0] = 6;
	edges_vertices[21][1] = 7;
	edges_vertices[22][0] = 7;
	edges_vertices[22][1] = 8;
	edges_vertices[23][0] = 8;
	edges_vertices[23][1] = 9;
	edges_vertices[24][0] = 9;
	edges_vertices[24][1] = 10;
	edges_vertices[25][0] = 6;
	edges_vertices[25][1] = 11;
	edges_vertices[26][0] = 7;
	edges_vertices[26][1] = 11;
	edges_vertices[27][0] = 8;
	edges_vertices[27][1] = 11;
	edges_vertices[28][0] = 9;
	edges_vertices[28][1] = 11;
	edges_vertices[29][0] = 10;
	edges_vertices[29][1] = 11;
	double vertices_edges[12][5];
	double face_vertices[20][3];
	face_vertices[0][0] = 0;
	face_vertices[0][1] = 1;
	face_vertices[0][2] = 2;
	face_vertices[1][0] = 0;
	face_vertices[1][1] = 2;
	face_vertices[1][2] = 3;
	face_vertices[2][0] = 0;
	face_vertices[2][1] = 3;
	face_vertices[2][2] = 4;
	face_vertices[3][0] = 0;
	face_vertices[3][1] = 4;
	face_vertices[3][2] = 5;
	face_vertices[4][0] = 0;
	face_vertices[4][1] = 5;
	face_vertices[4][2] = 1;
	face_vertices[5][0] = 1;
	face_vertices[5][1] = 10;
	face_vertices[5][2] = 6;
	face_vertices[6][0] = 6;
	face_vertices[6][1] = 2;
	face_vertices[6][2] = 1;
	face_vertices[7][0] = 2;
	face_vertices[7][1] = 6;
	face_vertices[7][2] = 7;
	face_vertices[8][0] = 7;
	face_vertices[8][1] = 3;
	face_vertices[8][2] = 2;
	face_vertices[9][0] = 3;
	face_vertices[9][1] = 7;
	face_vertices[9][2] = 8;
	face_vertices[10][0] = 8;
	face_vertices[10][1] = 3;
	face_vertices[10][2] = 4;
	face_vertices[11][0] = 4;
	face_vertices[11][1] = 9;
	face_vertices[11][2] = 8;
	face_vertices[12][0] = 9;
	face_vertices[12][1] = 4;
	face_vertices[12][2] = 5;
	face_vertices[13][0] = 5;
	face_vertices[13][1] = 10;
	face_vertices[13][2] = 9;
	face_vertices[14][0] = 10;
	face_vertices[14][1] = 1;
	face_vertices[14][2] = 5;
	face_vertices[15][0] = 11;
	face_vertices[15][1] = 6;
	face_vertices[15][2] = 10;
	face_vertices[16][0] = 11;
	face_vertices[16][1] = 7;
	face_vertices[16][2] = 6;
	face_vertices[17][0] = 11;
	face_vertices[17][1] = 7;
	face_vertices[17][2] = 8;
	face_vertices[18][0] = 11;
	face_vertices[18][1] = 8;
	face_vertices[18][2] = 9;
	face_vertices[19][0] = 11;
	face_vertices[19][1] = 9;
	face_vertices[19][2] = 10;
	double face_edges[20][3];
	face_edges[0][0] = 0;
	face_edges[0][1] = 5;
	face_edges[0][2] = 1;
	face_edges[1][0] = 1;
	face_edges[1][1] = 6;
	face_edges[1][2] = 2;
	face_edges[2][0] = 2;
	face_edges[2][1] = 7;
	face_edges[2][2] = 3;
	face_edges[3][0] = 3;
	face_edges[3][1] = 8;
	face_edges[3][2] = 4;
	face_edges[4][0] = 4;
	face_edges[4][1] = 9;
	face_edges[4][2] = 0;
	face_edges[5][0] = 19;
	face_edges[5][1] = 20;
	face_edges[5][2] = 10;
	face_edges[6][0] = 11;
	face_edges[6][1] = 5;
	face_edges[6][2] = 10;
	face_edges[7][0] = 11;
	face_edges[7][1] = 21;
	face_edges[7][2] = 12;
	face_edges[8][0] = 13;
	face_edges[8][1] = 6;
	face_edges[8][2] = 12;
	face_edges[9][0] = 13;
	face_edges[9][1] = 22;
	face_edges[9][2] = 14;
	face_edges[10][0] = 15;
	face_edges[10][1] = 7;
	face_edges[10][2] = 14;
	face_edges[11][0] = 15;
	face_edges[11][1] = 23;
	face_edges[11][2] = 16;
	face_edges[12][0] = 17;
	face_edges[12][1] = 8;
	face_edges[12][2] = 16;
	face_edges[13][0] = 17;
	face_edges[13][1] = 24;
	face_edges[13][2] = 18;
	face_edges[14][0] = 19;
	face_edges[14][1] = 9;
	face_edges[14][2] = 18;
	face_edges[15][0] = 25;
	face_edges[15][1] = 20;
	face_edges[15][2] = 29;
	face_edges[16][0] = 26;
	face_edges[16][1] = 21;
	face_edges[16][2] = 25;
	face_edges[17][0] = 27;
	face_edges[17][1] = 22;
	face_edges[17][2] = 26;
	face_edges[18][0] = 28;
	face_edges[18][1] = 23;
	face_edges[18][2] = 27;
	face_edges[19][0] = 29;
	face_edges[19][1] = 24;
	face_edges[19][2] = 28;
	double scale_height = 8e3;
	double atmos_height = scale_height*log(1+NUMBER_OF_LAYERS);
	for (int i = 0;i<=NUMBER_OF_SCALARS-1;++i)
	{
		double sigma, lat_edge_1, lon_edge_1, lat_edge_2, lon_edge_2, lat_point, lon_point, point_frac;
		int horizontal_index, edge_index,face_index;
		int layer_index = i/NUMBER_OF_SCALARS_H;
		sigma = 0.5*(scale_height/atmos_height)*(log((1+NUMBER_OF_LAYERS)/(layer_index+1))
		+log((1+NUMBER_OF_LAYERS)/(layer_index+2)));
		grid.z_scalar[i] = atmos_height*sigma;
		horizontal_index = i-(i/NUMBER_OF_SCALARS_H)*NUMBER_OF_SCALARS_H;
		if(horizontal_index < NUMBER_OF_PENTAGONS)
		{
			grid.latitude_scalar[i] = latitude_ico[horizontal_index];
			grid.longitude_scalar[i] = longitude_ico[horizontal_index];
			grid.adjacent_vector_indices_h[i][0] = NUMBER_OF_SCALARS_H + 
			layer_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_DUAL_VECTORS_H) + 5*horizontal_index;
			grid.adjacent_vector_indices_h[i][1] = NUMBER_OF_SCALARS_H + 
			layer_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_DUAL_VECTORS_H) + 5*horizontal_index+1;
			grid.adjacent_vector_indices_h[i][2] = NUMBER_OF_SCALARS_H + 
			layer_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_DUAL_VECTORS_H) + 5*horizontal_index+2;
			grid.adjacent_vector_indices_h[i][3] = NUMBER_OF_SCALARS_H + 
			layer_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_DUAL_VECTORS_H) + 5*horizontal_index+3;
			grid.adjacent_vector_indices_h[i][4] = NUMBER_OF_SCALARS_H + 
			layer_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_DUAL_VECTORS_H) + 5*horizontal_index+4;
			grid.adjacent_signs_h[i][0] = 1;
			grid.adjacent_signs_h[i][1] = 1;
			grid.adjacent_signs_h[i][2] = 1;
			grid.adjacent_signs_h[i][3] = 1;
			grid.adjacent_signs_h[i][4] = 1;
			grid.adjacent_signs_h[i][5] = 1;
			grid.vorticity_indices[i][0] = 0;
			grid.vorticity_indices[i][1] = 0;
			grid.vorticity_indices[i][2] = 0;
			grid.vorticity_indices[i][3] = 0;
			grid.vorticity_indices[i][4] = 0;
			grid.vorticity_indices[i][5] = 0;
			grid.vorticity_signs[i][0] = 1;
			grid.vorticity_signs[i][1] = 1;
			grid.vorticity_signs[i][2] = 1;
			grid.vorticity_signs[i][3] = 1;
			grid.vorticity_signs[i][4] = 1;
			grid.vorticity_signs[i][5] = 1;
		}
		else if(horizontal_index < NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES)
		{
			edge_index = (horizontal_index - NUMBER_OF_PENTAGONS)/POINTS_PER_EDGE;
			point_frac = (horizontal_index - NUMBER_OF_PENTAGONS -edge_index*POINTS_PER_EDGE)/POINTS_PER_EDGE;
			int vertic_index_1, vertic_index_2;
			vertic_index_1 = edges_vertices[edge_index][0];
			vertic_index_2 = edges_vertices[edge_index][1];
			find_geodetic(latitude_ico[vertic_index_1],longitude_ico[vertic_index_1],
			latitude_ico[vertic_index_2],longitude_ico[vertic_index_2],
			point_frac,&lat_point,&lon_point);
			grid.latitude_scalar[i] = lat_point;
			grid.longitude_scalar[i] = lon_point;
			grid.adjacent_vector_indices_h[i][0] = 0;
			grid.adjacent_vector_indices_h[i][1] = 0;
			grid.adjacent_vector_indices_h[i][2] = 0;
			grid.adjacent_vector_indices_h[i][3] = 0;
			grid.adjacent_vector_indices_h[i][4] = 0;
			grid.adjacent_signs_h[i][0] = 1;
			grid.adjacent_signs_h[i][1] = 1;
			grid.adjacent_signs_h[i][2] = 1;
			grid.adjacent_signs_h[i][3] = 1;
			grid.adjacent_signs_h[i][4] = 1;
			grid.adjacent_signs_h[i][5] = 1;
			grid.vorticity_indices[i][0] = 0;
			grid.vorticity_indices[i][1] = 0;
			grid.vorticity_indices[i][2] = 0;
			grid.vorticity_indices[i][3] = 0;
			grid.vorticity_indices[i][4] = 0;
			grid.vorticity_indices[i][5] = 0;
			grid.vorticity_signs[i][0] = 1;
			grid.vorticity_signs[i][1] = 1;
			grid.vorticity_signs[i][2] = 1;
			grid.vorticity_signs[i][3] = 1;
			grid.vorticity_signs[i][4] = 1;
			grid.vorticity_signs[i][5] = 1;
		}
		else
		{
			int inner_index = i - (NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES);
			face_index = (inner_index/POINTS_PER_INNER_FACES);
			int on_face_index = inner_index - face_index*POINTS_PER_INNER_FACES;
			int check, edge_1, edge_2, edge_index_1, edge_index_2;
			check = 1;
			edge_1 = face_edges[face_index][1];
			edge_2 = face_edges[face_index][2];
			int index_border, index_check, j, coord_1, coord_2, reverse_1, reverse_2;
			index_check = 0;
			coord_2 = -1;
			j = 0;
			while (check == 1)
			{
				++coord_2;
				index_border = pow(2,RES_ID) - 2 - j;
				if(j==0)
					--index_border;
				index_check = index_check + index_border;
				if(on_face_index <= index_check && on_face_index > index_check - index_border)
				{
					coord_1 = on_face_index - (index_check - index_border) - 1;
					check = 0;
				} 
				++j;
			}
			int index_1, index_2;
			reverse_1 = 1;
			reverse_2 = 1;
			if(face_vertices[face_index][1] == edges_vertices[edge_1][0])
				reverse_1 = 0;
			if(face_vertices[face_index][2] == edges_vertices[edge_2][0])
				reverse_2 = 0;
			if(reverse_1 == 0)
				index_1 = edge_1*POINTS_PER_EDGE + coord_2;
			else
				index_1 = edge_1*POINTS_PER_EDGE + POINTS_PER_EDGE - 1 - coord_2;
			if(reverse_2 == 0)
				index_2 = edge_2*POINTS_PER_EDGE + POINTS_PER_EDGE - 1 - coord_2;
			else
				index_2 = edge_2*POINTS_PER_EDGE + coord_2;
			double lat_1, lon_1, lat_2, lon_2, rel_on_line, lat_res, lon_res;
			lat_1 = grid.latitude_scalar[index_1];
			lon_1 = grid.longitude_scalar[index_1];
			lat_2 = grid.latitude_scalar[index_2];
			lon_2 = grid.longitude_scalar[index_2];
			rel_on_line = (index_border + 1 - coord_1)/(index_border + 2);
			find_geodetic(lat_1, lon_1, lat_2, lon_2, rel_on_line, &lat_res, &lon_res);
			grid.latitude_scalar[i] = lat_res;
			grid.longitude_scalar[i] = lon_res;
			grid.adjacent_vector_indices_h[i][0] = 0;
			grid.adjacent_vector_indices_h[i][1] = 0;
			grid.adjacent_vector_indices_h[i][2] = 0;
			grid.adjacent_vector_indices_h[i][3] = 0;
			grid.adjacent_vector_indices_h[i][4] = 0;
			grid.adjacent_vector_indices_h[i][5] = 0;
			grid.adjacent_signs_h[i][0] = 1;
			grid.adjacent_signs_h[i][1] = 1;
			grid.adjacent_signs_h[i][2] = 1;
			grid.adjacent_signs_h[i][3] = 1;
			grid.adjacent_signs_h[i][4] = 1;
			grid.adjacent_signs_h[i][5] = 1;
			grid.vorticity_indices[i][0] = 0;
			grid.vorticity_indices[i][1] = 0;
			grid.vorticity_indices[i][2] = 0;
			grid.vorticity_indices[i][3] = 0;
			grid.vorticity_indices[i][4] = 0;
			grid.vorticity_indices[i][5] = 0;
			grid.vorticity_signs[i][0] = 1;
			grid.vorticity_signs[i][1] = 1;
			grid.vorticity_signs[i][2] = 1;
			grid.vorticity_signs[i][3] = 1;
			grid.vorticity_signs[i][4] = 1;
			grid.vorticity_signs[i][5] = 1;
		}
		grid.adjacent_vector_index_upper[i] = i + layer_index*NUMBER_OF_VECTORS_H;
		grid.adjacent_vector_index_lower[i] = i + layer_index*NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_H + NUMBER_OF_SCALARS_H;
		double base_area = grid.area[grid.adjacent_vector_index_lower[i]];
		double radius_2 = semimajor + grid.z_vector[grid.adjacent_vector_index_upper[i]];
		double radius_1 = semimajor + grid.z_vector[grid.adjacent_vector_index_lower[i]];
		grid.volume[i] = (1/3)*(base_area/pow(radius_1,2))*(pow(radius_2,3)-pow(radius_1,3));
	}
	double sigma_1, sigma_2;
	for (int i=0;i<NUMBER_OF_VECTORS;++i)
	{
		int vert_index = i/(NUMBER_OF_SCALARS_H+NUMBER_OF_DUAL_VECTORS_H);
		int floor_index = vert_index*(NUMBER_OF_SCALARS_H+NUMBER_OF_DUAL_VECTORS_H);
		int h_index = i - floor_index;
		int on_edge_index;
		if(h_index >= NUMBER_OF_SCALARS_H)
		{
			if(h_index < NUMBER_OF_EDGES*(POINTS_PER_EDGE+1))
			{
				int edge_index = i/(POINTS_PER_EDGE+1);
				on_edge_index = i - edge_index*(POINTS_PER_EDGE+1);
				if(on_edge_index == 0)
				{
					grid.from_indices[i] = vertices_edges[edge_index][0];
					grid.to_indices[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE;
				}
				else if(on_edge_index == POINTS_PER_EDGE)
				{
					grid.from_indices[i] = vertices_edges[edge_index][1];
					grid.to_indices[i] = NUMBER_OF_PENTAGONS + (edge_index+1)*POINTS_PER_EDGE - 1;
				}
				else
				{
					grid.to_indices[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index;
					grid.from_indices[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index - 1;
				}	
			}
			else
			{
				int face_index = (h_index - NUMBER_OF_EDGES*(POINTS_PER_EDGE+1))/(VECTOR_POINTS_PER_INNER_FACE);
				int edge_1 = face_edges[face_index][1];
				int edge_2 = face_edges[face_index][2];
				int on_face_index =  NUMBER_OF_EDGES*(POINTS_PER_EDGE+1) + face_index*VECTOR_POINTS_PER_INNER_FACE;
				int triangle_on_face_index = on_face_index/3;
				int coord_1, coord_2;
				int check = 1;
				int index_border;
				int index_check = 0;
				coord_2 = -1;
				int j = 0;
				int reverse_1 = 1;
				int reverse_2 = 1;
				if(face_vertices[face_index][2] == edges_vertices[edge_2][0])
					reverse_2 = 0;
				if(face_vertices[face_index][1] == edges_vertices[edge_1][0])
					reverse_2 = 0;
				while (check == 1)
				{
					++coord_2;
					index_border = pow(2,RES_ID) - 1 - j;
					if(j==0)
						--index_border;
					index_check = index_check + index_border;
					if(on_face_index <= index_check && on_face_index > index_check - index_border)
					{
						coord_1 = triangle_on_face_index - (index_check - index_border) - 1;
						check = 0;
					} 
					++j;
				}
				int small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
				if(small_triangle_edge_index == 0)
				{
					if(coord_2 == 0)
					{
						grid.from_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES
						+ triangle_on_face_index - (index_border + 1);
						if(coord_1==0)
						{
							if(reverse_2 == 0)
							{
								grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][2] - 1
								+ POINTS_PER_EDGE - coord_2;
							}
							else
							{
								grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][1]
								+ coord_2;
							}
						}
						else
						{
							grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES
							+face_index*POINTS_PER_INNER_FACES + triangle_on_face_index - 1;
						}
					}
					else
					{
						grid.from_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES
						+ triangle_on_face_index - (index_border + 1);
						if(coord_1 == 0)
						{
							if(reverse_2 == 0)
							{
								grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][2]
								+ POINTS_PER_EDGE - 1 - coord_2;
							}
							else
							{
								grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][2]
								+ coord_2;
							}
						}
						else
						{
							grid.to_indices[i] = grid.from_indices[i] + index_border - 1;
						}
					}
				}
				if(small_triangle_edge_index == 1)
				{
					grid.to_indices[i] = 0;
					grid.from_indices[i] = 0;
				}
				if(small_triangle_edge_index == 2)
				{
					if(coord_1 == index_border - 1)
					{
						if(reverse_1 == 0)
						{
							grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][1]
							+ coord_2;
						}
						else
						{
							grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][1] + 
							POINTS_PER_EDGE - 1 - coord_2;
						}
					}
					else
					{
						grid.to_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES 
						+ triangle_on_face_index - coord_2;
					}
					if(coord_1 == 0)
					{
						if(reverse_2 == 0)
						{
							grid.from_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][2]
							+ coord_2;
						}
						else
						{
							grid.from_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][2]
							+ coord_2;
						}
					}
					else
					{
						grid.from_indices[i] = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES 
						+ triangle_on_face_index - coord_2 - 1;
					}
				}
			}
			double lat_point, lon_point;
			find_midpoint(grid.latitude_scalar[grid.from_indices[i] ],grid.longitude_scalar[grid.from_indices[i]],
			grid.latitude_scalar[grid.to_indices[i]],grid.longitude_scalar[grid.to_indices[i]],&lat_point,&lon_point);
			grid.latitude_vector[i] = lat_point;
			grid.longitude_vector[i] = lon_point;
			if(i < NUMBER_OF_SCALARS_H)
				sigma_1 = 1;
			else
				sigma_1 = (scale_height/atmos_height)*log((1+NUMBER_OF_LAYERS)/(vert_index+1));
			sigma_2 = (scale_height/atmos_height)*log((1+NUMBER_OF_LAYERS)/(vert_index+2));
			grid.gravity[i] = 0;
			grid.z_vector[i] = atmos_height*0.5*(sigma_1+sigma_2);
			grid.parallel_distance[i] = 0;
			grid.normal_distance[i] = calculate_distance_h(grid.latitude_vector[grid.to_indices[i]],
			grid.longitude_vector[grid.to_indices[i]],grid.latitude_vector[grid.from_indices[i]],
			grid.longitude_vector[grid.from_indices[i]],semimajor+grid.z_vector[i]);
			grid.area[i] = calculate_vertical_face(grid.parallel_distance[i],semimajor+atmos_height*sigma_2,
			semimajor+atmos_height*sigma_1);
		}
		else
		{
			grid.to_indices[i] = h_index + (vert_index-1)*NUMBER_OF_SCALARS_H;
			grid.from_indices[i] = h_index + vert_index*NUMBER_OF_SCALARS_H;
			if(i>=NUMBER_OF_VECTORS_H+NUMBER_OF_SCALARS_H)
			{
				grid.latitude_vector[i] = grid.latitude_scalar[grid.to_indices[i]];
				grid.longitude_vector[i] = grid.longitude_scalar[grid.to_indices[i]];
			}
			else
			{
				grid.latitude_vector[i] = grid.latitude_scalar[grid.from_indices[i]];
				grid.longitude_vector[i] = grid.longitude_scalar[grid.from_indices[i]];
			}
			if(i < NUMBER_OF_SCALARS_H)
				sigma_1 = 1;
			else
				sigma_1  = grid.z_scalar[grid.to_indices[i]]/atmos_height;
			sigma_2 = grid.z_scalar[grid.from_indices[i]]/atmos_height;
			int sigma = (scale_height/atmos_height)*log((1+NUMBER_OF_LAYERS)/(vert_index+1));
			grid.z_vector[i] = atmos_height*sigma;
			if(h_index < NUMBER_OF_PENTAGONS)
				grid.area[i] = 5*4*M_PI*pow(semimajor+grid.z_vector[i],2)/(3*NUMBER_OF_TRIANGLES);
			else
				grid.area[i] = 2*4*M_PI*pow(semimajor+grid.z_vector[i],2)/NUMBER_OF_TRIANGLES;
			grid.normal_distance[i] = fabs(atmos_height*sigma_1-atmos_height*sigma_2);
			grid.gravity[i] = -9.8;
		}
	}
	for (int i=0;i<NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H;++i)
	{
		grid.h_curl_indices[i][0] = 1;
		grid.h_curl_indices[i][1] = 1;
		grid.h_curl_indices[i][2] = 1;
		grid.h_curl_indices[i][3] = 1;
		grid.h_curl_signs[i][0] = 1;
		grid.h_curl_signs[i][1] = 1;
		grid.h_curl_signs[i][2] = 1;
		grid.h_curl_signs[i][3] = 1;
		grid.vector_product_sign[i] = 1;
		grid.recov_hor_par_dual_index[i][0] = 0;
		grid.recov_hor_par_dual_index[i][1] = 0;
		grid.recov_hor_par_dual_weight[i][0] = 0;
		grid.recov_hor_par_dual_weight[i][1] = 0;
		grid.recov_hor_par_pri_index[i][0] = 0;
		grid.recov_hor_par_pri_index[i][1] = 0;
		grid.recov_hor_par_pri_index[i][2] = 0;
		grid.recov_hor_par_pri_index[i][3] = 0;
		grid.recov_hor_par_pri_index[i][4] = 0;
		grid.recov_hor_par_pri_index[i][5] = 0;
		grid.recov_hor_par_pri_index[i][6] = 0;
		grid.recov_hor_par_pri_index[i][7] = 0;
		grid.recov_hor_par_pri_index[i][8] = 0;
		grid.recov_hor_par_pri_index[i][9] = 0;
		grid.recov_hor_par_pri_index[i][10] = 0;
		grid.recov_hor_par_pri_weight[i][0] = 0;
		grid.recov_hor_par_pri_weight[i][1] = 0;
		grid.recov_hor_par_pri_weight[i][2] = 0;
		grid.recov_hor_par_pri_weight[i][3] = 0;
		grid.recov_hor_par_pri_weight[i][4] = 0;
		grid.recov_hor_par_pri_weight[i][5] = 0;
		grid.recov_hor_par_pri_weight[i][6] = 0;
		grid.recov_hor_par_pri_weight[i][7] = 0;
		grid.recov_hor_par_pri_weight[i][8] = 0;
		grid.recov_hor_par_pri_weight[i][9] = 0;
		grid.recov_hor_par_pri_weight[i][10] = 0;
		grid.recov_hor_ver_dual_index[i][0] = 0;
		grid.recov_hor_ver_dual_index[i][1] = 0;
		grid.recov_hor_ver_dual_weight[i][0] = 0;
		grid.recov_hor_ver_dual_weight[i][1] = 0;
		grid.recov_hor_ver_pri_index[i][0] = 0;
		grid.recov_hor_ver_pri_index[i][1] = 0;
		grid.recov_hor_ver_pri_index[i][2] = 0;
		grid.recov_hor_ver_pri_index[i][3] = 0;
		grid.recov_hor_ver_pri_weight[i][0] = 0;
		grid.recov_hor_ver_pri_weight[i][1] = 0;
		grid.recov_hor_ver_pri_weight[i][2] = 0;
		grid.recov_hor_ver_pri_weight[i][3] = 0;
	}
	for (int i=0;i<(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H;++i)
	{
		grid.recov_ver_1_dual_index[i][0] = 0;
		grid.recov_ver_1_dual_index[i][1] = 0;
		grid.recov_ver_1_dual_index[i][2] = 0;
		grid.recov_ver_1_dual_index[i][3] = 0;
		grid.recov_ver_1_dual_index[i][4] = 0;
		grid.recov_ver_1_dual_index[i][5] = 0;
		grid.recov_ver_1_dual_weight[i][0] = 0;
		grid.recov_ver_1_dual_weight[i][1] = 0;
		grid.recov_ver_1_dual_weight[i][2] = 0;
		grid.recov_ver_1_dual_weight[i][3] = 0;
		grid.recov_ver_1_dual_weight[i][4] = 0;
		grid.recov_ver_1_dual_weight[i][5] = 0;
		grid.recov_ver_1_pri_index[i][0] = 0;
		grid.recov_ver_1_pri_index[i][1] = 0;
		grid.recov_ver_1_pri_index[i][2] = 0;
		grid.recov_ver_1_pri_index[i][3] = 0;
		grid.recov_ver_1_pri_index[i][4] = 0;
		grid.recov_ver_1_pri_index[i][5] = 0;
		grid.recov_ver_1_pri_weight[i][0] = 0;
		grid.recov_ver_1_pri_weight[i][1] = 0;
		grid.recov_ver_1_pri_weight[i][2] = 0;
		grid.recov_ver_1_pri_weight[i][3] = 0;
		grid.recov_ver_1_pri_weight[i][4] = 0;
		grid.recov_ver_1_pri_weight[i][5] = 0;
		grid.recov_ver_2_dual_index[i][0] = 0;
		grid.recov_ver_2_dual_index[i][1] = 0;
		grid.recov_ver_2_dual_index[i][2] = 0;
		grid.recov_ver_2_dual_index[i][3] = 0;
		grid.recov_ver_2_dual_index[i][4] = 0;
		grid.recov_ver_2_dual_index[i][5] = 0;
		grid.recov_ver_2_dual_weight[i][0] = 0;
		grid.recov_ver_2_dual_weight[i][1] = 0;
		grid.recov_ver_2_dual_weight[i][2] = 0;
		grid.recov_ver_2_dual_weight[i][3] = 0;
		grid.recov_ver_2_dual_weight[i][4] = 0;
		grid.recov_ver_2_dual_weight[i][5] = 0;
		grid.recov_ver_2_pri_index[i][0] = 0;
		grid.recov_ver_2_pri_index[i][1] = 0;
		grid.recov_ver_2_pri_index[i][2] = 0;
		grid.recov_ver_2_pri_index[i][3] = 0;
		grid.recov_ver_2_pri_index[i][4] = 0;
		grid.recov_ver_2_pri_index[i][5] = 0;
		grid.recov_ver_2_pri_weight[i][0] = 0;
		grid.recov_ver_2_pri_weight[i][1] = 0;
		grid.recov_ver_2_pri_weight[i][2] = 0;
		grid.recov_ver_2_pri_weight[i][3] = 0;
		grid.recov_ver_2_pri_weight[i][4] = 0;
		grid.recov_ver_2_pri_weight[i][5] = 0;
	}
	for (int i=0;i<NUMBER_OF_DUAL_SCALARS;++i)
	{
		int layer_index = i/NUMBER_OF_DUAL_SCALARS_H;
		int h_index = i - layer_index*NUMBER_OF_DUAL_SCALARS_H;
		int face_index = h_index/TRIANGLES_PER_FACE;
		int edge_1 = face_edges[face_index][1];
		int edge_2 = face_edges[face_index][2];
		int on_face_index =  h_index - face_index*TRIANGLES_PER_FACE;
		int index_1, index_2, index_3;
		int coord_1, coord_2;
		int check = 1;
		int index_border;
		int index_check = 0;
		coord_2 = -1;
		int j = 0;
		int reverse_1 = 1;
		int reverse_2 = 1;
		if(face_vertices[face_index][2] == edges_vertices[edge_2][0])
			reverse_2 = 0;
		if(face_vertices[face_index][1] == edges_vertices[edge_1][0])
			reverse_2 = 0;
		int points_right = 0;
		while (check == 1)
		{
			if(points_right == 0)
				points_right = 1;
			if(points_right == 1)
				points_right = 0;
			++coord_2;
			index_border = pow(2,RES_ID) - j;
			if(j==0)
				--index_border;
			index_check = index_check + index_border;
			if(on_face_index <= index_check && on_face_index > index_check - index_border)
			{
				coord_1 = on_face_index - (index_check - index_border) - 1;
				check = 0;
			} 
			++j;
		}
		int column_index = coord_2/2;
		if(points_right == 1)
		{
			if(column_index == 0)
			{
				if(coord_1 == 0)
					index_1 = face_vertices[face_index][0];
				else
					index_1 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*face_edges[face_index][0] + coord_1 - 1;
			}
			else
			{
				if(coord_1 == 0)
				{
					if(reverse_2 == 0)
					{
						;
					}
					else
					{
						;
					}
				}
				else
				{
					index_1 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*TRIANGLES_PER_FACE
					+ 0;
				}
			}
		}
		else
		{
			;
		}
		double lat_res, lon_res;
		find_center(grid.latitude_scalar[index_1],grid.longitude_scalar[index_1],
		grid.latitude_scalar[index_2],grid.longitude_scalar[index_2],
		grid.latitude_scalar[index_3],grid.longitude_scalar[index_3],&lat_res,&lon_res);
		dualgrid.latitude_scalar[i] = lat_res;
		dualgrid.longitude_scalar[i] = lon_res;
		dualgrid.z_scalar[i]
		= atmos_height*0.5*(scale_height/atmos_height)*(log((1+NUMBER_OF_LAYERS)/(layer_index+1))
		+log((1+NUMBER_OF_LAYERS)/(layer_index+2)));
		dualgrid.adjacent_vector_indices_h[i][0] = 0;
		dualgrid.adjacent_vector_indices_h[i][1] = 0;
		dualgrid.adjacent_vector_indices_h[i][2] = 0;
		dualgrid.adjacent_vector_index_lower[i] = 0;
		dualgrid.adjacent_vector_index_upper[i] = 0;
		dualgrid.adjacent_signs_h[i][0] = 1;
		dualgrid.adjacent_signs_h[i][1] = 1;
		dualgrid.adjacent_signs_h[i][2] = 1;
	}
	for (int i=0;i<NUMBER_OF_DUAL_VECTORS;++i)
	{
		dualgrid.longitude_vector[i] = 0;
		dualgrid.latitude_vector[i] = 0;
		int layer_index = i/NUMBER_OF_DUAL_VECTORS_H;
		dualgrid.z_vector[i] = atmos_height*(scale_height/atmos_height)*log((1+NUMBER_OF_LAYERS)/(layer_index+1));
		dualgrid.direction[i] = 0;
		double z_1, z_2;
		if(i == 0)
			z_1 = atmos_height;
		else
			z_1 = atmos_height*0.5*(scale_height/atmos_height)*(log((1+NUMBER_OF_LAYERS)/(i))+log((1+NUMBER_OF_LAYERS)/(i+1)));
		z_2 = atmos_height*0.5*(scale_height/atmos_height)*(log((1+NUMBER_OF_LAYERS)/(i+1))+log((1+NUMBER_OF_LAYERS)/(i+2)));
		if (dualgrid.is_vertical[i] > 0)
			dualgrid.f_vec[i] = 2*omega*sin(dualgrid.latitude_vector[i]);
		else
			dualgrid.f_vec[i] = 2*omega*cos(dualgrid.latitude_vector[i])*cos(dualgrid.direction[i]);
	}
	for (int i=0;i<NUMBER_OF_TRIANGLES*NUMBER_OF_LAYERS;++i)
	{
		dualgrid.vorticity_indices[i][0] = 0;
		dualgrid.vorticity_indices[i][1] = 0;
		dualgrid.vorticity_indices[i][2] = 0;
		dualgrid.vorticity_signs[i][0] = 1;
		dualgrid.vorticity_signs[i][1] = 1;
		dualgrid.vorticity_signs[i][2] = 1;
	}
	for (int i=0;i<NUMBER_OF_DUAL_VECTORS_H*(NUMBER_OF_LAYERS+1);++i)
	{
		dualgrid.h_curl_indices[i][0] = 0;
		dualgrid.h_curl_indices[i][1] = 0;
		dualgrid.h_curl_indices[i][2] = 0;
		dualgrid.h_curl_indices[i][3] = 0;
		dualgrid.h_curl_signs[i][0] = 1;
		dualgrid.h_curl_signs[i][1] = 1;
		dualgrid.h_curl_signs[i][2] = 1;
		dualgrid.h_curl_signs[i][3] = 1;
	}
	find_min_dist();
}

void find_min_dist(void)
{
	for (int i = 0; i<NUMBER_OF_VECTORS; ++i)
	{
		double dist = grid.normal_distance[i];
		if (dist < min_dist)
			min_dist = dist;
	}
}

void find_geos(double x,double y,double z,double *lat_out,double *lon_out)
{
	*lat_out = asin(z/sqrt(x*x+y*y+z*z));
	*lon_out = atan2(y,x);
	if(lon_out < 0)
		*lon_out = *lon_out + 2*M_PI;
}

double calculate_distance_h(double latitude_a,double longitude_a,double latitude_b,double longitude_b,double radius)
{
	double dist = 2*radius*asin(sqrt(1/2-(1/2)*(cos(latitude_a)*cos(latitude_b)*cos(longitude_b-longitude_a)
	+sin(latitude_a)*sin(latitude_b))));
	return dist;
}
 
double calculate_vertical_face(double base_distance,double r_1,double r_2)
{
	double area;
	area = base_distance*(0.5*pow(r_2,2)/r_1 - r_1/2);
	return area;
}
 
void find_midpoint(double lat_1_in,double lon_1_in,double lat_2_in,double lon_2_in,double *lat_out,double *lon_out)
{
	double x_1, y_1, z_1, x_2, y_2, z_2, x, y, z;
	x_1 = cos(lat_1_in)*cos(lon_1_in);
	y_1 = cos(lat_1_in)*sin(lon_1_in);
	z_1 = sin(lat_1_in);
	x_2 = cos(lat_2_in)*cos(lon_2_in);
	y_2 = cos(lat_2_in)*sin(lon_2_in);
	z_2 = sin(lat_2_in);
	x = x_1 + x_2;
	y = y_1 + y_2;
	z = z_1 + z_2;
	find_geos(x,y,z,lat_out,lon_out);
}

void find_geodetic(double lat_1_in,double lon_1_in,double lat_2_in,double lon_2_in,
double parameter,double *lat_out,double *lon_out)
{
	double tau_dash = 0.5 + (0.5+0.5*(cos(lat_1_in)*cos(lat_2_in)*cos(lon_1_in-lon_2_in)+sin(lat_1_in)*sin(lat_2_in)))*
	tan(2*M_PI*parameter-asin(sqrt(0.5-0.5*(cos(lat_1_in)*cos(lat_2_in)*cos(lon_1_in-lon_2_in)+sin(lat_1_in)*sin(lat_2_in)))));
	*lat_out = asin(tau_dash*sin(lat_2_in)+(1-tau_dash)*sin(lat_1_in));
	*lon_out = atan2(tau_dash*cos(lat_2_in)*sin(lon_2_in)+(1-tau_dash)*cos(lat_1_in)*sin(lon_1_in),
	tau_dash*cos(lat_2_in)*cos(lon_2_in)+(1-tau_dash)*cos(lat_1_in)*cos(lon_1_in));
	if(*lon_out < 0)
		*lon_out = *lon_out + 2*M_PI;
}

double find_geodetic_direction(double lat_1_in,double lon_1_in,double lat_2_in,double lon_2_in,
double parameter)
{
	double rel_vec[3], local_i[3], local_j[3];
	rel_vec[0] = cos(lat_2_in)*cos(lon_2_in) - cos(lat_1_in)*cos(lon_1_in);
	rel_vec[1] = cos(lat_2_in)*sin(lon_2_in) - cos(lat_1_in)*sin(lon_1_in);
	rel_vec[2] = sin(lat_2_in) - sin(lat_1_in);
	double lat, lon;
	lat = 0;
	lon = 0;
	find_geodetic(lat_1_in,lon_1_in,lat_2_in,lon_2_in,parameter,&lat,&lon);
	local_i[0] = -sin(lon);
	local_i[1] = cos(lon);
	local_i[2] = 0;
	local_j[0] = -sin(lat)*cos(lon);
	local_j[1] = -sin(lat)*sin(lon);
	local_j[2] = cos(lat);
	double x_comp = scalar_product_simple(local_i,rel_vec);
	double y_comp = scalar_product_simple(local_j,rel_vec);
	double direction = atan2(x_comp,y_comp);
	if(direction < 0)
		direction = direction + 2*M_PI;
	return direction;
}

void find_center(double lat_1_in,double lon_1_in,double lat_2_in,double lon_2_in,double lat_3_in,double lon_3_in,
double *lat_out,double *lon_out)
{
	double epsilon = semimajor*pow(0.5,RES_ID+10);
	double lat_1,lon_1,lat_2,lon_2,lat_3,lon_3;
	lat_1 = lat_1_in;
	lon_1 = lon_1_in;
	lat_2 = lat_2_in;
	lon_2 = lon_2_in;
	lat_3 = lat_3_in;
	lon_3 = lon_3_in;
	double dist_1, dist_2, dist_3;
	dist_1 = calculate_distance_h(lat_1_in,lon_1_in,lat_2_in,lon_2_in,semimajor);
	dist_2 = calculate_distance_h(lat_1_in,lon_1_in,lat_3_in,lon_3_in,semimajor);
	dist_3 = calculate_distance_h(lat_2_in,lon_2_in,lat_3_in,lon_3_in,semimajor);
	double max_dist = dist_1;
	if(dist_2 > dist_1)
		max_dist = dist_2;
	if(dist_3 > max_dist)
		max_dist = dist_3;
	double lat_int_1, lon_int_1, lat_int_2, lon_int_2, lat_int_3, lon_int_3;
	while(max_dist >= epsilon)
	{
		find_midpoint(lat_1,lon_1,lat_2,lon_2,&lat_int_1,&lon_int_1);
		find_midpoint(lat_1,lon_1,lat_3,lon_3,&lat_int_2,&lon_int_2);
		find_midpoint(lat_2,lon_2,lat_3,lon_3,&lat_int_3,&lon_int_3);
		dist_1 = calculate_distance_h(lat_1_in,lon_1_in,lat_2_in,lon_2_in,semimajor);
		dist_2 = calculate_distance_h(lat_1_in,lon_1_in,lat_3_in,lon_3_in,semimajor);
		dist_3 = calculate_distance_h(lat_2_in,lon_2_in,lat_3_in,lon_3_in,semimajor);
		double max_dist = dist_1;
		if(dist_2 > dist_1)
			max_dist = dist_2;
		if(dist_3 > max_dist)
			max_dist = dist_3;
	}
	*lat_out = lat_1;
	*lon_out = lon_1;
}
