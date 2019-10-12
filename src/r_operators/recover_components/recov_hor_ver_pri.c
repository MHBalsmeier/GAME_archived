#include <math.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../r_operators.h"


void recov_hor_ver_pri(Vector_field in_field, double out_field[])
{
	extern Grid grid;
	for(int i = 0; i<NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H; ++i)
	{
		out_field[i] = 
		grid.recov_hor_ver_pri_weight[i][0]*in_field[grid.recov_hor_ver_pri_index[i][0]] + 
		grid.recov_hor_ver_pri_weight[i][1]*in_field[grid.recov_hor_ver_pri_index[i][1]] + 
		grid.recov_hor_ver_pri_weight[i][2]*in_field[grid.recov_hor_ver_pri_index[i][2]] + 
		grid.recov_hor_ver_pri_weight[i][3]*in_field[grid.recov_hor_ver_pri_index[i][3]];
	}
}
