#include "../enum_and_typedefs.h"

int dissipation(Vector_field velocity, Scalar_field density, Vector_field friction_acc, Scalar_field heating, Grid *grid)
{
	for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
	{
		friction_acc[i] = 0;
	}
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	{
		heating[i] = 0;
	}
	return 0;
}
