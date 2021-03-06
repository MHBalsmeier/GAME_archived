/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"

int vector_times_vector(Vector_field in_field_0, Vector_field in_field_1, Vector_field out_field)
{
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	out_field[i] = in_field_0[i]*in_field_1[i];
	}
    return 0;
}
