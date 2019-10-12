#include <math.h>
#include "../enum_and_typedefs.h"
#include "io.h"

State interpolation_t(State state_m1, State state_0, float t_m1, float t_0, float t_write)
{
	double weight_m1, weight_0;
	weight_0 = (t_write-t_m1)/(t_0-t_m1);
	weight_m1 = 1 - weight_0;
	State state_write;
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		state_write.density[i] = weight_m1*state_m1.density[i] + weight_0*state_0.density[i];
		state_write.pot_temp[i] = weight_m1*state_m1.pot_temp[i] + weight_0*state_0.pot_temp[i];
	}
	for (int i = 0; i<=NUMBER_OF_VECTORS-1; ++i)
	{
		state_write.wind[i] = weight_m1*state_m1.wind[i] + weight_0*state_0.wind[i];

	}
	return state_write;
}
