/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains the (ocean surface) wave component of GAME.
*/

#include "../enum_and_typedefs.h"

int update_waves(Waves *waves, Diagnostics *diag, Grid *grid)
{
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		if (grid -> is_land[i] == 0)
		{
			waves -> swh[i] =  0.25*pow(diag -> e_kin[NO_OF_SCALARS - NO_OF_SCALARS_H + i], 0.5);
		}
	}
	return 0;
}
