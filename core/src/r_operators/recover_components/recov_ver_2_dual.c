#include "../../enum_and_typedefs.h"

void recov_ver_2_dual(Dual_vector_field in_field, double out_field[])
{
    extern Grid grid;
    for(int i = 0; i < NUMBER_OF_V_VECTORS; ++i)
    {
        for (int j = 0; j < 6; j++)
            out_field[i] = out_field[i] + grid.recov_ver_2_dual_weight[6*i + j]*in_field[grid.recov_ver_2_dual_index[6*i + j]]; 
    }
}
