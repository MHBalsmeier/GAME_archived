#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "r_operators.h"


void coriolis(Vector_field wind, Vector_field acc)
{
	extern Grid grid;
	extern Dualgrid dualgrid;
	vector_product(wind,dualgrid.f_vec,acc);
}
