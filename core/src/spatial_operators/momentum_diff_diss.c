/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#include <stdlib.h>
#include <stdio.h>

int momentum_diff_diss(State *state, Diagnostics *diagnostics, Irreversible_quantities *irrev, Config_info *config_info, Grid *grid, Dualgrid *dualgrid, double delta_t)
{
	int layer_index, h_index;
	// Firstly the diffusion.
	// Evaluating necessary differential operators.
	divv_h(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
	if (config_info -> momentum_diff_v == 1)
	{
		add_vertical_divv(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
    }
	if (config_info -> momentum_diff_v == 0)
	{
    	grad_hor(diagnostics -> velocity_gas_divv, irrev -> velocity_grad_div, grid);
	}
	if (config_info -> momentum_diff_v == 1)
	{
    	grad(diagnostics -> velocity_gas_divv, irrev -> velocity_grad_div, grid);
	}
    curl_of_vorticity_m(diagnostics -> rel_vort, diagnostics -> curl_of_vorticity_m, grid, dualgrid, config_info);
    // Multiplying the values of the differential operators by the effective viscosities.
    // Calculating the effective viscosity of the divergence term.
    calc_divv_term_viscosity_eff(state, irrev -> viscosity_eff, grid, delta_t);
	scalar_times_vector(irrev -> viscosity_eff, irrev -> velocity_grad_div, irrev -> friction_acc, grid);
    // Calculating the effective viscosity of the curl term.
    calc_curl_term_viscosity_eff(state, irrev -> viscosity_eff, grid, delta_t);
	scalar_times_vector(irrev -> viscosity_eff, diagnostics -> curl_of_vorticity_m, diagnostics -> curl_of_vorticity_m, grid);
	// adding the curl term to the divergence term
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		irrev -> friction_acc[i] += diagnostics -> curl_of_vorticity_m[i];
	}
	
	// at these very coarse resolutions, a divergence damping must be added to control grid-scale noise (switched out currently)
	if (RES_ID <= -1)
	{
		// the order of the divergence damping
		int div_damp_order = 4;
		// checking if the divergence damping order is even
		if (fmod(div_damp_order, 2) != 0)
		{
			printf("div_damp_order must be even.");
			printf("Aborting.\n");
			exit(1);
		}
		divv_h(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
		grad(diagnostics -> velocity_gas_divv, irrev -> velocity_grad_div, grid);
    	for (int i = 0; i < div_damp_order/2 - 1; ++i)
    	{
			// diagnostics -> velocity_gas_divv is a misuse of name
			divv_h(irrev -> velocity_grad_div, diagnostics -> velocity_gas_divv, grid);
			// irrev -> velocity_grad_div is a misuse of name here, it is actually the divegence damping
    		grad(diagnostics -> velocity_gas_divv, irrev -> velocity_grad_div, grid);
		}
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_H_VECTORS; ++i)
		{
			layer_index = i/NO_OF_VECTORS_H;
			h_index = i - layer_index*NO_OF_VECTORS_H;
			irrev -> friction_acc[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index]
			// a sign has to be taken into account here
			+= -pow(-1, div_damp_order/2)*config_info -> div_damp_coeff*irrev -> velocity_grad_div[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
		}
	}
	
	// simplified dissipation
	inner_product(state -> velocity_gas, irrev -> friction_acc, irrev -> heating_diss, grid);
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irrev -> heating_diss[i] = -density_gas(state, i)*irrev -> heating_diss[i];
	}
	return 0;
}

int curl_of_vorticity_m(Curl_field vorticity, Vector_field out_field, Grid *grid, Dualgrid *dualgrid, Config_info *config_info)
{
	// Calculates the negative curl of the vorticity.
	int layer_index, h_index, no_of_edges, i, j;
	#pragma omp parallel for private(layer_index, h_index, no_of_edges, i, j)
	for (i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		out_field[i] = 0;
		// The vertical case.
		if (h_index < NO_OF_SCALARS_H)
		{
			no_of_edges = 6;
			if (h_index < NO_OF_PENTAGONS)
			{
				no_of_edges = 5;
			}
			for (j = 0; j < no_of_edges; ++j)
			{
				out_field[i] += grid -> adjacent_signs_h[6*h_index + j]
				*dualgrid -> normal_distance[layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
				*vorticity[layer_index*2*NO_OF_VECTORS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
			}
			// Dividing by the area. The additional minus sign accounts for the fact the we want the negative curl of the curl.
			out_field[i] = -config_info -> momentum_diff_v*out_field[i]/grid -> area[i];
		}
		// The horizontal case.
		// Remember: (curl(zeta))*e_x = dzeta_z/dy - dzeta_y/dz = (dz*dzeta_z - dy*dzeta_y)/(dy*dz) = (dz*dzeta_z - dy*dzeta_y)/area (Stokes' Theorem, which is used here)
		else
		{
			// horizontal difference of vertical vorticity (dzeta_z*dz)
			// An averaging over three rhombi must be done.
			for (j = 0; j < 3; ++j)
			{
				out_field[i] +=
				// This prefactor accounts for the fact that we average over three rhombi.
				1.0/3*(
				// vertical length at the to_index_dual point
				dualgrid -> normal_distance[NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + dualgrid -> to_index[h_index - NO_OF_SCALARS_H]]
				// vorticity at the to_index_dual point
				*vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + dualgrid -> adjacent_vector_indices_h[3*dualgrid -> to_index[h_index - NO_OF_SCALARS_H] + j]]
				// vertical length at the from_index_dual point
				- dualgrid -> normal_distance[NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + dualgrid -> from_index[h_index - NO_OF_SCALARS_H]]
				// vorticity at the from_index_dual point
				*vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + dualgrid -> adjacent_vector_indices_h[3*dualgrid -> from_index[h_index - NO_OF_SCALARS_H] + j]]);
			}
			if (config_info -> momentum_diff_v == 1)
			{
				// vertical difference of horizontal vorticity (-dzeta_y*dy)
				out_field[i] +=
				// substracting the upper zeta_y value
				-dualgrid -> normal_distance[layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + h_index - NO_OF_SCALARS_H]
				*vorticity[layer_index*2*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
				out_field[i] +=
				// adding the lower zeta_y value
				dualgrid -> normal_distance[(layer_index + 1)*NO_OF_DUAL_VECTORS_PER_LAYER + h_index - NO_OF_SCALARS_H]
				*vorticity[(layer_index + 1)*2*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
			}
			// Dividing by the area. The additional minus sign accounts for the fact the we want the negative curl of the curl.
			out_field[i] = -out_field[i]/grid -> area[i];
		}
	}
	return 0;
}










