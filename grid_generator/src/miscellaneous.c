/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"
#include <math.h>

int set_f_vec(double latitude_vector[], double direction_dual[], double latitude_vector_dual[], double f_vec[])
{
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H; ++i)
    {
        if (i >= NUMBER_OF_DUAL_VECTORS_H)
        	f_vec[i] = 2*OMEGA*sin(latitude_vector[i - NUMBER_OF_DUAL_VECTORS_H]);
        else
        	f_vec[i] = 2*OMEGA*cos(latitude_vector_dual[i])*sin(direction_dual[i]);
    }
 	return 0;   
}
