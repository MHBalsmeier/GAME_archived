#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"
#include <stdio.h>

int temperature_diagnostics(Scalar_field density_entropy, Scalar_field density, Scalar_field temperature)
{
    double exner_pressure, pot_temp;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	pot_temp = exp(density_entropy[i]/(C_P*density[i]) - K_B*N_A/(M_D*C_P)*log(1/P_0*K_B*K_B*pow(M_D/N_A*exp(5.0/3)/(M_PI*H_BAR*H_BAR), 1.5)));
        exner_pressure = pow(R_D*density[i]*pot_temp/P_0, R_D/C_V);
        temperature[i] = exner_pressure*pot_temp;
    }
    return 0;
}
