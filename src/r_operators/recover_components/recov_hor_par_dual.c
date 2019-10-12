#include <math.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../r_operators.h"


void recov_hor_par_dual(Dual_vector_field in_field, double out_field[])
{
	extern Grid grid;
	for(int i = 0; i<NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H; ++i)
	{
		out_field[i] = 
		grid.recov_hor_par_dual_weight[i][0]*in_field[grid.recov_hor_par_dual_index[i][0]] + 
		grid.recov_hor_par_dual_weight[i][1]*in_field[grid.recov_hor_par_dual_index[i][1]];
	}
}
