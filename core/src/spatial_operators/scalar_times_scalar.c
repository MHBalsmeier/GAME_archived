/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"

int scalar_times_scalar(Scalar_field in_field_0, Scalar_field in_field_1, Scalar_field out_field)
{
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	out_field[i] = in_field_0[i]*in_field_1[i];
	}
    return 0;
}
