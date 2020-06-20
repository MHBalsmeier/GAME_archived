#include "../enum_and_typedefs.h"
#include <stdio.h>

int kinetic_energy(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        out_field[i] = 0;
        for (int j = 0; j < 14; ++j)
        	out_field[i] += grid -> e_kin_weights[14*i + j]*pow(in_field[grid -> e_kin_indices[14*i + j]], 2);
    }
    return 0;
}
