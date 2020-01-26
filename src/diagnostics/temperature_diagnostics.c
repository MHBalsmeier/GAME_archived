#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"

void temperature_diagnostics(Scalar_field pot_temp, Scalar_field pressure, Scalar_field temperature)
{
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
		temperature[i] = pot_temp[i]*pow(pressure[i]/P_0 ,R_D/C_P);
}
