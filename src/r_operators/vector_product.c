#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"

void vector_product (Vector_field a_field, Dual_vector_field b_field, Vector_field out_field)
{
	extern Grid grid;
	extern Dualgrid dualgrid;
	double vector_1, vector_2, vector_3, vector_4;
	double hor_par_pri[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H];
	double hor_par_dual[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H];
	double hor_ver_pri[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H];
	double hor_ver_dual[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H];
	double ver_1_pri[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H];
	double ver_1_dual[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H];
	double ver_2_pri[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H];
	double ver_2_dual[(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H];
	recov_hor_par_dual(b_field,hor_par_dual);
	recov_hor_par_pri(a_field,hor_par_pri);
	recov_hor_ver_dual(b_field,hor_ver_dual);
	recov_hor_ver_pri(a_field,hor_ver_pri);
	recov_ver_1_dual(b_field,ver_1_dual);
	recov_ver_1_pri(a_field,ver_1_pri);
	recov_ver_2_dual(b_field,ver_2_dual);
	recov_ver_2_pri(a_field,ver_2_pri);
	for (int i = 0; i<=NUMBER_OF_VECTORS-1; ++i)
	{
		double term_1, term_2;
		if(grid.is_vertical[i] == 0)
		{
			int sign;
			recov_hor_par_dual(b_field,hor_par_dual);
			recov_hor_par_pri(a_field,hor_par_pri);
			recov_hor_ver_dual(b_field,hor_ver_dual);
			recov_hor_ver_pri(a_field,hor_ver_pri);
			vector_1 = hor_par_pri[grid.vertical_horizontal_index[i]];
			vector_2 = hor_par_dual[grid.vertical_horizontal_index[i]];
			vector_3 = hor_ver_pri[grid.vertical_horizontal_index[i]];
			vector_4 = hor_ver_dual[grid.vertical_horizontal_index[i]];
			sign = grid.vector_product_sign[grid.vertical_horizontal_index[i]];
			term_1 = sign*vector_1*vector_4;
			term_2 = -sign*vector_2*vector_3;
		}
		else
		{
			recov_ver_1_pri(a_field,ver_1_pri);
			recov_ver_1_dual(b_field,ver_1_dual);
			recov_ver_2_pri(a_field,ver_2_pri);
			recov_ver_2_dual(b_field,ver_2_dual);
			vector_1 = ver_1_pri[grid.vertical_horizontal_index[i]];
			vector_2 = ver_1_dual[grid.vertical_horizontal_index[i]];
			vector_3 = ver_2_pri[grid.vertical_horizontal_index[i]];
			vector_4 = ver_2_dual[grid.vertical_horizontal_index[i]];
			term_1 = vector_1*vector_2;
			term_2 = -vector_3*vector_4;
		}
		out_field[i] = term_1, term_2;
	}
}
