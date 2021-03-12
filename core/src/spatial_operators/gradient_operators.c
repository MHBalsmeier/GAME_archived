/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include "spatial_operators.h"
#include <stdio.h>

int grad_hor_cov(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	// calculates the horizontal covariant gradient
    int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            out_field[i]
            = (in_field[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]
            - in_field[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H])
            /grid -> normal_distance[i];
        }
    }
    return 0;
}

int grad_vert_cov(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	// calculates the vertical covariant gradient
    int layer_index, h_index, lower_index, upper_index;
    // loop over the inner grid points
	#pragma omp parallel for private(layer_index, h_index, lower_index, upper_index)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_SCALARS_H)
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            out_field[i] = (in_field[upper_index] - in_field[lower_index])/grid -> normal_distance[i];
        }
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        out_field[i] =
        out_field[i + NO_OF_VECTORS_PER_LAYER]
        + (out_field[i + NO_OF_VECTORS_PER_LAYER] - out_field[i + 2*NO_OF_VECTORS_PER_LAYER])/
        (grid -> z_vector[i + NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[i + 2*NO_OF_VECTORS_PER_LAYER])
        *(grid -> z_vector[i] - grid -> z_vector[i + NO_OF_VECTORS_PER_LAYER]);
    }
    // linear extrapolation to the surface
    #pragma omp parallel for
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        out_field[i] =
        out_field[i - NO_OF_VECTORS_PER_LAYER]
        + (out_field[i - 2*NO_OF_VECTORS_PER_LAYER] - out_field[i - NO_OF_VECTORS_PER_LAYER])/
        (grid -> z_vector[i - 2*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[i - NO_OF_VECTORS_PER_LAYER])
        *(grid -> z_vector[i] - grid -> z_vector[i - NO_OF_VECTORS_PER_LAYER]);
    }
    return 0;
}

int grad_oro_corr(Vector_field con_to_cov_gradient_field, Grid *grid)
{
    // correction for terrain
    int layer_index, h_index;
    double vertical_gradient;
	#pragma omp parallel for private(layer_index, h_index, vertical_gradient)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H && layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
        {
        	remap_verpri2horpri_vector(con_to_cov_gradient_field, layer_index, h_index - NO_OF_SCALARS_H, &vertical_gradient, grid);
            con_to_cov_gradient_field[i] += -grid -> slope[i]*vertical_gradient;
        }
    }
    return 0;
}

int grad_cov(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	grad_hor_cov(in_field, out_field, grid);
	grad_vert_cov(in_field, out_field, grid);
    return 0;
}

int grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	grad_cov(in_field, out_field, grid);
	grad_oro_corr(out_field, grid);
    return 0;
}

int grad_hor(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	grad_cov(in_field, out_field, grid);
	grad_oro_corr(out_field, grid);
    int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_SCALARS_H)
        {
            out_field[i] = 0;
        }
    }
    return 0;
}




