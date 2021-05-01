/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
In this file, the gravity potential gets computed.
*/

#include "include.h"
#include "enum.h"
#include "geos95.h"

int set_gravity_potential(double z_scalar[], double gravity_potential[], double GRAVITY_MEAN_SFC_ABS)
{
	/*
	This function computes the gravity potential.
	*/
	
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	gravity_potential[i] = -GRAVITY_MEAN_SFC_ABS*(RADIUS*RADIUS/(RADIUS + z_scalar[i]) - RADIUS);
    }
	return 0;
}
