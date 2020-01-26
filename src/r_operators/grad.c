#include "../enum_and_typedefs.h"

void grad(Scalar_field in_field, Vector_field out_field)
{
	extern Grid grid;
    for (int i = NUMBER_OF_SCALARS_H; i < NUMBER_OF_VECTORS - NUMBER_OF_SCALARS_H; i++)
        out_field[i] = (in_field[grid.to_indices[i]] - in_field[grid.from_indices[i]])/grid.normal_distance[i];
    for (int i = 0; i < NUMBER_OF_SCALARS_H; i++)
        out_field[i] = out_field[i + NUMBER_OF_VECTORS_H + NUMBER_OF_SCALARS_H];
    for (int i = NUMBER_OF_VECTORS - NUMBER_OF_SCALARS_H; i < NUMBER_OF_VECTORS; i++)
        out_field[i] = out_field[i - (NUMBER_OF_VECTORS_H + NUMBER_OF_SCALARS_H)];
}
