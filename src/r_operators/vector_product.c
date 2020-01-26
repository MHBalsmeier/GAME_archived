#include "../enum_and_typedefs.h"
#include "r_operators.h"

void vector_product(Vector_field a_field, Dual_vector_field b_field, Vector_field out_field)
{
	extern Grid grid;
	extern Dualgrid dualgrid;
	double vector_1,  vector_2,  vector_3,  vector_4, term_1, term_2;
    int sign;
	double hor_par_pri[NUMBER_OF_H_VECTORS];
	double hor_par_dual[NUMBER_OF_H_VECTORS];
	double hor_ver_pri[NUMBER_OF_H_VECTORS];
	double hor_ver_dual[NUMBER_OF_H_VECTORS];
	double ver_1_pri[NUMBER_OF_V_VECTORS];
	double ver_1_dual[NUMBER_OF_V_VECTORS];
	double ver_2_pri[NUMBER_OF_V_VECTORS];
	double ver_2_dual[NUMBER_OF_V_VECTORS];
	recov_hor_par_dual(b_field, hor_par_dual);
	recov_hor_par_pri(a_field, hor_par_pri);
	recov_hor_ver_dual(b_field, hor_ver_dual);
	recov_hor_ver_pri(a_field, hor_ver_pri);
	recov_ver_1_dual(b_field, ver_1_dual);
	recov_ver_1_pri(a_field, ver_1_pri);
	recov_ver_2_dual(b_field, ver_2_dual);
	recov_ver_2_pri(a_field, ver_2_pri);
    int vertical_horizontal_index, vert_index, floor_index, h_index;
	for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
	{
        vert_index = i/(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
		floor_index = vert_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
		h_index = i - floor_index;
        vertical_horizontal_index = 0;
		if(h_index >= NUMBER_OF_SCALARS_H)
		{
            
			recov_hor_par_dual(b_field, hor_par_dual);
			recov_hor_par_pri(a_field, hor_par_pri);
			recov_hor_ver_dual(b_field, hor_ver_dual);
			recov_hor_ver_pri(a_field, hor_ver_pri);
			vector_1 = hor_par_pri[vertical_horizontal_index];
			vector_2 = hor_par_dual[vertical_horizontal_index];
			vector_3 = hor_ver_pri[vertical_horizontal_index];
			vector_4 = hor_ver_dual[vertical_horizontal_index];
			sign = grid.vector_product_sign[vertical_horizontal_index];
			term_1 = sign*vector_1*vector_4;
			term_2 = -sign*vector_2*vector_3;
		}
		else
		{
			recov_ver_1_pri(a_field, ver_1_pri);
			recov_ver_1_dual(b_field, ver_1_dual);
			recov_ver_2_pri(a_field, ver_2_pri);
			recov_ver_2_dual(b_field, ver_2_dual);
			vector_1 = ver_1_pri[vertical_horizontal_index];
			vector_2 = ver_1_dual[vertical_horizontal_index];
			vector_3 = ver_2_pri[vertical_horizontal_index];
			vector_4 = ver_2_dual[vertical_horizontal_index];
			term_1 = vector_1*vector_2;
			term_2 = -vector_3*vector_4;
		}
		out_field[i] = term_1 + term_2;
	}
}
