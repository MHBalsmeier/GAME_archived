#include <math.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../r_operators.h"


void recov_hor_par_pri(Vector_field in_field, double out_field[])
{
	extern Grid grid;
	for(int i = 0; i<NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H; ++i)
	{
		out_field[i] = 
		grid.recov_hor_par_pri_weight[i][0]*in_field[grid.recov_hor_par_pri_index[i][0]] + 
		grid.recov_hor_par_pri_weight[i][1]*in_field[grid.recov_hor_par_pri_index[i][1]] + 
		grid.recov_hor_par_pri_weight[i][2]*in_field[grid.recov_hor_par_pri_index[i][2]] +
		grid.recov_hor_par_pri_weight[i][3]*in_field[grid.recov_hor_par_pri_index[i][3]] + 
		grid.recov_hor_par_pri_weight[i][4]*in_field[grid.recov_hor_par_pri_index[i][4]] + 
		grid.recov_hor_par_pri_weight[i][5]*in_field[grid.recov_hor_par_pri_index[i][5]] + 
		grid.recov_hor_par_pri_weight[i][6]*in_field[grid.recov_hor_par_pri_index[i][6]] + 
		grid.recov_hor_par_pri_weight[i][7]*in_field[grid.recov_hor_par_pri_index[i][7]] + 
		grid.recov_hor_par_pri_weight[i][8]*in_field[grid.recov_hor_par_pri_index[i][8]] +
		grid.recov_hor_par_pri_weight[i][9]*in_field[grid.recov_hor_par_pri_index[i][9]] + 
		grid.recov_hor_par_pri_weight[i][10]*in_field[grid.recov_hor_par_pri_index[i][10]];
	}
}
