/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"

int momentum_diff_diss(Vector_field velocity, Scalar_field density, Vector_field friction_acc, Scalar_field heating, Config_info *config_info, Grid *grid)
{
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		friction_acc[i] = 0;
	}
	// The heating is a power density, as such it has the SI unit J/(m^3s).
	inner_product(velocity, friction_acc, heating, grid);
	scalar_times_scalar(density, heating, heating);
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		heating[i] = -C_D_V*heating[i];
	}
	return 0;
}
