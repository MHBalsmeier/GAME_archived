#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"
#include <stdio.h>

int exner_pressure_diagnostics(Scalar_field density_pot_temp, Scalar_field exner_pressure)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        exner_pressure[i] = pow(R_D*density_pot_temp[i]/P_0, R_D/C_V);
    return 0;
}
