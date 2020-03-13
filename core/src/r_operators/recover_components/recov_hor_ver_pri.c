#include "../../enum_and_typedefs.h"

void recov_hor_ver_pri(Vector_field in_field, double out_field[], Grid *grid)
{
    for(int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        out_field[i] = 0;
        for (int j = 0; j < 4; ++j)
            out_field[i] = out_field[i] + grid -> recov_hor_ver_pri_weight[4*i + j]*in_field[grid -> recov_hor_ver_pri_index[4*i + j]];
    }
}
