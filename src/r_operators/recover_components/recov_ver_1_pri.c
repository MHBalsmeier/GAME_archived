#include <math.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../r_operators.h"


void recov_ver_1_pri(Vector_field in_field, double out_field[])
{
	extern Grid grid;
	for(int i = 0; i<(NUMBER_OF_LAYERS+1)*NUMBER_OF_SCALARS_H; ++i)
	{
		out_field[i] = 
		grid.recov_ver_1_pri_weight[i][0]*in_field[grid.recov_ver_1_pri_index[i][0]] + 
		grid.recov_ver_1_pri_weight[i][1]*in_field[grid.recov_ver_1_pri_index[i][1]] + 
		grid.recov_ver_1_pri_weight[i][2]*in_field[grid.recov_ver_1_pri_index[i][2]] + 
		grid.recov_ver_1_pri_weight[i][3]*in_field[grid.recov_ver_1_pri_index[i][3]] + 
		grid.recov_ver_1_pri_weight[i][4]*in_field[grid.recov_ver_1_pri_index[i][4]] + 
		grid.recov_ver_1_pri_weight[i][5]*in_field[grid.recov_ver_1_pri_index[i][5]];
	}
}
