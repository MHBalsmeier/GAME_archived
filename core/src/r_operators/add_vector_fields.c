#include "../enum_and_typedefs.h"

int add_vector_fields(Vector_field a_field, Vector_field b_field, Vector_field out_field)
{
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
        out_field[i] = a_field[i] + b_field[i];
    return 0;
}
