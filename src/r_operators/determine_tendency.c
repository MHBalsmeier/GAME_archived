#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"

State tendency(State current_state)
{
	State state_tendency;
	extern Grid grid;
	Scalar_field adv_density, adv_pot_temp, v_divergence;
	adv_s(current_state.wind, current_state.density,adv_density);
	adv_s(current_state.wind, current_state.pot_temp,adv_pot_temp);
	divergence(current_state.wind, v_divergence);
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		state_tendency.density[i] = -adv_density[i] - current_state.density[i]*v_divergence[i];
		state_tendency.pot_temp[i] = -adv_pot_temp[i];
	}
	Vector_field adv_momentum, pressure_gradient, coriolis_tendency;
	Scalar_field pressure;
	pressure_diagnostics(current_state.pot_temp,current_state.density,pressure);
	grad(pressure,pressure_gradient);
	adv_v(current_state.wind,adv_momentum);
	coriolis(current_state.wind,coriolis_tendency);
	for (int i = 0; i<=NUMBER_OF_VECTORS-1; ++i)
	{
		state_tendency.wind[i] = -adv_momentum[i] - (1/current_state.density[i])*pressure_gradient[i]
		+ coriolis_tendency[i] + grid.gravity[i];
		if(i <=NUMBER_OF_VECTORS-1-NUMBER_OF_VECTORS_H)
		{
			state_tendency.wind[i] = 0;
		}
	}
	return state_tendency;
}
