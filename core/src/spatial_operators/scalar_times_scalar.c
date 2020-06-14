#include "../enum_and_typedefs.h"

int scalar_times_scalar(Scalar_field in_field_0, Scalar_field in_field_1, Scalar_field out_field)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    	out_field[i] = in_field_0[i]*in_field_1[i];
    return 0;
}
