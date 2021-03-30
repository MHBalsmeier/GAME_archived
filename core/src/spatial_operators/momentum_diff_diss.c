/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int hori_momentum_diffusion(State *state, Diagnostics *diagnostics, Irreversible_quantities *irrev, Config_info *config_info, Grid *grid, Dualgrid *dualgrid, double delta_t)
{
	// Evaluating necessary differential operators.
	divv_h(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
	grad_hor(diagnostics -> velocity_gas_divv, irrev -> velocity_grad_div, grid);
    curl_of_vorticity(diagnostics -> rel_vort, diagnostics -> curl_of_vorticity, grid, dualgrid, config_info);
	// adding the curl term to the divergence term
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		irrev -> friction_acc[i] = irrev -> velocity_grad_div[i] - diagnostics -> curl_of_vorticity[i];
	}
    // Calculating the effective horizontal kinematic viscosity (Eddy viscosity).
	hori_viscosity_eff(state, irrev -> viscosity_eff, grid, diagnostics, config_info, delta_t);
	// multiplying by the viscosity
	vector_times_vector(irrev -> viscosity_eff, irrev -> friction_acc, irrev -> friction_acc);
	return 0;
}

int vert_momentum_diffusion(State *state, Irreversible_quantities *irrev, Grid *grid)
{
	// some parameters
	double bndr_lr_height = 1e3; // boundary layer height
	double bndr_lr_visc_max = 1.0/86400; // maximum friction coefficient in the boundary layer
	double e_folding_height = 0.5*bndr_lr_height;
	int layer_index, h_index, vector_index;
	double z_agl;
	#pragma omp parallel for private(layer_index, h_index, vector_index, z_agl)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
		// height above ground level
		z_agl = grid -> z_vector[vector_index]
		- 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[h_index]]);
		// adding the boundary layer friction
		if (z_agl < bndr_lr_height)
		{
			irrev -> friction_acc[vector_index]
			+= -bndr_lr_visc_max*(exp(-z_agl/e_folding_height) - exp(-bndr_lr_height/e_folding_height))
			/(1 - exp(-bndr_lr_height/e_folding_height))
			*state -> velocity_gas[vector_index];
		}
	}
	return 0;
}

int curl_of_vorticity(Curl_field vorticity, Vector_field out_field, Grid *grid, Dualgrid *dualgrid, Config_info *config_info)
{
	// Calculates the negative curl of the vertical vorticity.
	int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		out_field[i] = 0;
		// Remember: (curl(zeta))*e_x = dzeta_z/dy - dzeta_y/dz = (dz*dzeta_z - dy*dzeta_y)/(dy*dz) = (dz*dzeta_z - dy*dzeta_y)/area (Stokes' Theorem, which is used here)
		if (h_index >= NO_OF_SCALARS_H)
		{
			// horizontal difference of vertical vorticity (dzeta_z*dz)
			// An averaging over three rhombi must be done.
			for (int j = 0; j < 3; ++j)
			{
				out_field[i] = out_field[i]
				// This prefactor accounts for the fact that we average over three rhombi.
				+ 1.0/3*(
				// vertical length at the to_index_dual point
				dualgrid -> normal_distance[NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + dualgrid -> to_index[h_index - NO_OF_SCALARS_H]]
				// vorticity at the to_index_dual point
				*vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + dualgrid -> adjacent_vector_indices_h[3*dualgrid -> to_index[h_index - NO_OF_SCALARS_H] + j]]
				// vertical length at the from_index_dual point
				- dualgrid -> normal_distance[NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + dualgrid -> from_index[h_index - NO_OF_SCALARS_H]]
				// vorticity at the from_index_dual point
				*vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + dualgrid -> adjacent_vector_indices_h[3*dualgrid -> from_index[h_index - NO_OF_SCALARS_H] + j]]);
			}
		}
		// Dividing by the area.
		out_field[i] = out_field[i]/grid -> area[i];
	}
	return 0;
}

int simple_dissipation_rate(State *state, Irreversible_quantities *irrev, Grid *grid)
{
	// simplified dissipation
	inner_product(state -> velocity_gas, irrev -> friction_acc, irrev -> heating_diss, grid, 1);
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irrev -> heating_diss[i] = -density_gas(state, i)*irrev -> heating_diss[i];
	}
	return 0;
}








