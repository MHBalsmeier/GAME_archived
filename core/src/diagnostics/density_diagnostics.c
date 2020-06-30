/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <stdio.h>
#include <stdlib.h>
#include "../enum_and_typedefs.h"

double calc_micro_density(double density_macro, double condensates_density_sum)
{
	double result = density_macro/(1 - condensates_density_sum/RHO_WATER);
	if (result < -EPSILON_TRACERS/(1 - condensates_density_sum/RHO_WATER))
	{
		printf("Error: microscopic density negative.\n");
		exit(1);
	}
	return result;
}

double calc_condensates_density_sum(int layer_index, int h_index, Tracer_densities tracer_densities)
{
	double result = 0;
	for (int i = 0; i < NO_OF_CONDENSATED_TRACERS; ++i)
		result += tracer_densities[i*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + h_index];
	if (result < -NO_OF_CONDENSATED_TRACERS*EPSILON_TRACERS)
	{
		printf("Error: condensates_density_sum negative.\n");
		exit(1);
	}
	if (result >= RHO_WATER)
	{
		printf("Error: condensates_density_sum >= RHO_WATER.\n");
		exit(1);
	}
	return result;
}