#include "../enum_and_typedefs.h"

void divergence(Vector_field in_field, Scalar_field out_field)
{
	extern Grid grid;
    long layer_index, horizontal_index;
    double comp_h, comp_v;
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	{
        out_field[i] = 0;
        layer_index = i/NUMBER_OF_SCALARS_H;
		horizontal_index = i - layer_index*NUMBER_OF_SCALARS_H;
		if (horizontal_index < NUMBER_OF_PENTAGONS)
		{
            for (int j = 0; j < 5; j++)
                comp_h += in_field[grid.adjacent_vector_indices_h[6*i + j]]*grid.adjacent_signs_h[6*i + j]*grid.area[grid.adjacent_vector_indices_h[6*i + j]];
		}
		else
		{
            for (int j = 0; j < 6; j++)
                comp_h += in_field[grid.adjacent_vector_indices_h[6*i + j]]*grid.adjacent_signs_h[6*i + j]*grid.area[grid.adjacent_vector_indices_h[6*i + j]];
		}
        comp_v = in_field[grid.adjacent_vector_index_upper[i]]*grid.area[grid.adjacent_vector_index_upper[i]] - in_field[grid.adjacent_vector_index_lower[i]]*grid.area[grid.adjacent_vector_index_lower[i]];
		out_field[i] = (1/grid.volume[i])*(comp_h + comp_v);
	}
}
