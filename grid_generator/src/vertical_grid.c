/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int determine_z_scalar(double z_scalar[], double z_vertical_vector_pre[], double z_surface[], double z_oro_off, double TOA)
{
	int layer_index, h_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
		h_index = i - layer_index*NUMBER_OF_SCALARS_H;
		for (int j = 0; j < NUMBER_OF_LAYERS + 1; ++j)
		{
			if (j >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
				z_vertical_vector_pre[j] = z_surface[h_index] + (z_oro_off - z_surface[h_index])/NUMBER_OF_ORO_LAYERS*(NUMBER_OF_LAYERS - j);
			else
				z_vertical_vector_pre[j] = TOA - (TOA - z_oro_off)/(NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)*j;
		}
		z_scalar[i] = 0.5*(z_vertical_vector_pre[layer_index] + z_vertical_vector_pre[layer_index + 1]);
        if (z_scalar[i] <= 0)
		{
            printf("z_scalar contains a non-positive value.\n");
			exit(1);
		}
    }
    return 0;
}
