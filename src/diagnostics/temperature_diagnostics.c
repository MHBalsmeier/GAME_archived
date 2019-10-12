#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"

void temperature_diagnostics(Scalar_field pot_temp, Scalar_field pressure, Scalar_field temperature)
{
	extern double R_d, p_0, c_p;
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		temperature[i] = pot_temp[i]*pow(pressure[i]/p_0,R_d/c_p);
	}
}
