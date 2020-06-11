#include <stdio.h>
#include <stdlib.h>
#include "../enum_and_typedefs.h"

double calc_micro_density(double density_macro, double condensates_density_sum)
{
	double epsilon = 0.00001;
	double result = density_macro/(1 - condensates_density_sum/RHO_WATER);
	if (result < -epsilon)
	{
		printf("Error: microscopic density negative.\n");
		exit(1);
	}
	return result;
}

double calc_condensates_density_sum(int layer_index, int h_index, Add_comp_densities add_comp_densities)
{
	double result = 0;
	double epsilon = 0.00001;
	for (int i = 0; i < NUMBER_OF_COND_ADD_COMPS; ++i)
		result += add_comp_densities[i*NUMBER_OF_SCALARS + layer_index*NUMBER_OF_SCALARS_H + h_index];
	if (result < -epsilon)
	{
		printf("Error: condensates_density_sum negative.\n");
		exit(1);
	}
	return result;
}
