/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdio.h>

int scalar_times_vector(Scalar_field scalar_field, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
	/*
	This function multiplies the vector field vector_field by the scalar field scalar_field.
	*/
    int layer_index, h_index, lower_index, upper_index, i;
    double scalar_value;
	#pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, i, scalar_value)
    for (i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value
            = 0.5*(
            scalar_field[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]
            + scalar_field[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            scalar_value = 0.5*(
            scalar_field[upper_index]
            + scalar_field[lower_index]);
        }
    	out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value
        = scalar_field[i]
        + (scalar_field[i] - scalar_field[i + NO_OF_SCALARS_H])
        /(grid -> z_scalar[i] - grid -> z_scalar[i + NO_OF_SCALARS_H])
        *(grid -> z_vector[i] - grid -> z_scalar[i]);
    	out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the surface
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value
        = scalar_field[upper_index]
        + (scalar_field[upper_index - NO_OF_SCALARS_H] - scalar_field[upper_index])
        /(grid -> z_scalar[upper_index - NO_OF_SCALARS_H] - grid -> z_scalar[upper_index])
        *(grid -> z_vector[i] - grid -> z_scalar[upper_index]);
    	out_field[i] = scalar_value*vector_field[i];
    }
    return 0;
}

int advection_3rd_order(Scalar_field tracer, Vector_field mass_flux_density, Vector_field wind, Vector_field gradient, Vector_field out_field, Grid *grid, int rk_substep, double delta_t)
{
	/*
	This function calculates third order scalar tracer advection.
	*/
	grad_hor(tracer, gradient, grid);
    int layer_index, h_index, lower_index, upper_index, i;
    // determining the sign for the upstream correction
    int sign = -1 + 2*rk_substep;
    double tracer_value, tangential_wind_value, tangential_gradient_value;
	#pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, i, tracer_value, tangential_wind_value, tangential_gradient_value)
    for (i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        // 3rd order advection is only done in the horizontal
        if (h_index >= NO_OF_SCALARS_H)
        {
        	tangential_wind(wind, layer_index, h_index - NO_OF_SCALARS_H, &tangential_wind_value, grid);
        	tangential_wind(gradient, layer_index, h_index - NO_OF_SCALARS_H, &tangential_gradient_value, grid);
            tracer_value
            = 0.5*(
            tracer[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]
            + tracer[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
            // adding the upstream correction
            tracer_value += sign*0.5*delta_t*
            (gradient[i]*wind[i]
            + tangential_gradient_value*tangential_wind_value);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            tracer_value = 0.5*(
            tracer[upper_index]
            + tracer[lower_index]);
        }
    	out_field[i] = tracer_value*mass_flux_density[i];
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for private(tracer_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        tracer_value
        = tracer[i]
        + (tracer[i] - tracer[i + NO_OF_SCALARS_H])
        /(grid -> z_scalar[i] - grid -> z_scalar[i + NO_OF_SCALARS_H])
        *(grid -> z_vector[i] - grid -> z_scalar[i]);
    	out_field[i] = tracer_value*mass_flux_density[i];
    }
    // linear extrapolation to the surface
    #pragma omp parallel for private(layer_index, h_index, upper_index, tracer_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        tracer_value
        = tracer[upper_index]
        + (tracer[upper_index - NO_OF_SCALARS_H] - tracer[upper_index])
        /(grid -> z_scalar[upper_index - NO_OF_SCALARS_H] - grid -> z_scalar[upper_index])
        *(grid -> z_vector[i] - grid -> z_scalar[upper_index]);
    	out_field[i] = tracer_value*mass_flux_density[i];
    }
    return 0;
}

int scalar_times_vector_scalar_h_v(Scalar_field in_field_h, Scalar_field in_field_v, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index, i;
    double scalar_value;
    #pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, i, scalar_value)
    for (i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value
            = 0.5*(
            in_field_h[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]
            + in_field_h[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            scalar_value = 0.5*(
            in_field_v[upper_index]
            + in_field_v[lower_index]);
        }
        out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for private(scalar_value)
    for (i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value
        = in_field_v[i]
        + (in_field_v[i] - in_field_v[i + NO_OF_SCALARS_H])
        /(grid -> z_scalar[i] - grid -> z_scalar[i + NO_OF_SCALARS_H])
        *(grid -> z_vector[i] - grid -> z_scalar[i]);
        out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the surface
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value
        = in_field_v[upper_index]
        + (in_field_v[upper_index - NO_OF_SCALARS_H] - in_field_v[upper_index])
        /(grid -> z_scalar[upper_index - NO_OF_SCALARS_H] - grid -> z_scalar[upper_index])
        *(grid -> z_vector[i] - grid -> z_scalar[upper_index]);
        out_field[i] = scalar_value*vector_field[i];
    }
    return 0;
}





