#include "../enum_and_typedefs.h"
#include "r_operators.h"
#include <stdlib.h>
#include <stdio.h>

int vector_product(Vector_field a_field, Dual_vector_field b_field, Vector_field out_field, Grid *grid)
{
    double vector_1, vector_2, vector_3, vector_4, term_1, term_2;
    short sign;
    double *hor_par_pri = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *hor_par_dual = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *hor_ver_pri = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *hor_ver_dual = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *ver_1_pri = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *ver_1_dual = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *ver_2_pri = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    double *ver_2_dual = malloc(NUMBER_OF_H_VECTORS*sizeof(double));
    recov_hor_par_dual(b_field, hor_par_dual, grid);
    recov_hor_par_pri(a_field, hor_par_pri, grid);
    recov_hor_ver_dual(b_field, hor_ver_dual, grid);
    recov_hor_ver_pri(a_field, hor_ver_pri, grid);
    recov_ver_1_dual(b_field, ver_1_dual, grid);
    recov_ver_1_pri(a_field, ver_1_pri, grid);
    recov_ver_2_dual(b_field, ver_2_dual, grid);
    recov_ver_2_pri(a_field, ver_2_pri, grid);
    int vertical_horizontal_index, vert_index, floor_index, h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        vert_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - vert_index*NUMBER_OF_VECTORS_PER_LAYER;
        vertical_horizontal_index = 0;
        if(h_index >= NUMBER_OF_VECTORS_V)
        {
            recov_hor_par_dual(b_field, hor_par_dual, grid);
            recov_hor_par_pri(a_field, hor_par_pri, grid);
            recov_hor_ver_dual(b_field, hor_ver_dual, grid);
            recov_hor_ver_pri(a_field, hor_ver_pri, grid);
            vector_1 = hor_par_pri[vertical_horizontal_index];
            vector_2 = hor_par_dual[vertical_horizontal_index];
            vector_3 = hor_ver_pri[vertical_horizontal_index];
            vector_4 = hor_ver_dual[vertical_horizontal_index];
            sign = grid -> vector_product_sign[vertical_horizontal_index];
            term_1 = sign*vector_1*vector_4;
            term_2 = -sign*vector_2*vector_3;
        }
        else
        {
            recov_ver_1_pri(a_field, ver_1_pri, grid);
            recov_ver_1_dual(b_field, ver_1_dual, grid);
            recov_ver_2_pri(a_field, ver_2_pri, grid);
            recov_ver_2_dual(b_field, ver_2_dual, grid);
            vector_1 = ver_1_pri[vertical_horizontal_index];
            vector_2 = ver_1_dual[vertical_horizontal_index];
            vector_3 = ver_2_pri[vertical_horizontal_index];
            vector_4 = ver_2_dual[vertical_horizontal_index];
            term_1 = vector_1*vector_2;
            term_2 = -vector_3*vector_4;
        }
        out_field[i] = term_1 + term_2;
    }
    free(hor_par_pri);
    free(hor_par_dual);
    free(hor_ver_pri);
    free(hor_ver_dual);
    free(ver_1_pri);
    free(ver_1_dual);
    free(ver_2_pri);
    free(ver_2_dual);
    return 0;
}
