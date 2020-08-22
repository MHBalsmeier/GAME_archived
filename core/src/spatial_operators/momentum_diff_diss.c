/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"

int momentum_diff_diss(Vector_field velocity, Scalar_field density, Vector_field friction_acc, Scalar_field heating, Config_info *config_info, Grid *grid)
{
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		friction_acc[i] = 0;
	}
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		heating[i] = 0;
	}
	return 0;
}
