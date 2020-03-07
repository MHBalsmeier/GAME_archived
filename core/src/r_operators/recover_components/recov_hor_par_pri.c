#include "../../enum_and_typedefs.h"

void recov_hor_par_pri(Vector_field in_field, double out_field[])
{
    extern Grid grid;
    for( int i = 0; i < NUMBER_OF_H_VECTORS; ++i)
    {
        for (int j = 0; j < 2; j++)
            out_field[i] = out_field[i] + grid.recov_hor_par_pri_weight[2*i + j]*in_field[grid.recov_hor_par_pri_index[2*i + j]];
    }
}
