#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "time_stepping.h"
#include "../r_operators/r_operators.h"

State leapfrog(State state_m2, State state_m1, State tendency_m1, float delta_t)
{
	extern double rho_water;
	State state_0;
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		state_0.density[i] = state_m2.density[i] + 2*delta_t*tendency_m1.density[i];
		state_0.pot_temp[i] = state_m2.pot_temp[i] + 2*delta_t*tendency_m1.pot_temp[i];
	}
	for (int i = 0; i<=NUMBER_OF_VECTORS-1; ++i)
	{
		state_0.wind[i] = state_m2.wind[i] + 2*delta_t*tendency_m1.wind[i];
	}
	Scalar_field pressure;
	pressure_diagnostics(state_0.pot_temp,state_0.density,pressure);
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		state_0.pressure[i] = pressure[i];
	}
	return state_0;
}
