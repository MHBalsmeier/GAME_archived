#include "../../enum_and_typedefs.h"

void recov_ver_1_pri(Vector_field in_field, double out_field[], Grid *grid)
{
    for(int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        out_field[i] = 0;
        for (int j = 0; j < 6; ++j)
            out_field[i] = out_field[i] + grid -> recov_ver_1_pri_weight[6*i + j]*in_field[grid -> recov_ver_1_pri_index[6*i + j]]; 
    }
}
