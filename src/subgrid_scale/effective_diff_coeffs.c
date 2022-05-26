/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, diffusion coefficients, including Eddy viscosities, are computed.
*/

#include <math.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "subgrid_scale.h"

double tke2vert_diff_coeff(double);
double tke2hor_diff_coeff(double);

int hor_div_viscosity(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config *config)
{
	/*
	This function computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal divergent movements.
	*/
	
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irrev -> viscosity_div[i] = density_gas(state, i)*tke2hor_diff_coeff(irrev -> tke[i]);
	}
	return 0;
}

int hor_curl_viscosity_rhombi(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config *config)
{
	/*
	This function computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal curl movements on rhombi.
	*/
	double molecular_viscosity;
	
	int scalar_index_from, scalar_index_to, vector_index;
	#pragma omp parallel for private(molecular_viscosity, scalar_index_from, scalar_index_to, vector_index)
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
			
			// indices of the adjacent scalar grid points
			scalar_index_from = layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index];
			scalar_index_to = layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index];
			
			// preliminary result
			irrev -> viscosity_curl_rhombi[vector_index] = 0.5*(tke2hor_diff_coeff(irrev -> tke[scalar_index_from]) + tke2hor_diff_coeff(irrev -> tke[scalar_index_to]));
			
			// calculating and adding the molecular viscosity
			molecular_viscosity = 0.5*(irrev -> molecular_diffusion_coeff[scalar_index_from] + irrev -> molecular_diffusion_coeff[scalar_index_to]);
			irrev -> viscosity_curl_rhombi[vector_index] += molecular_viscosity;
			
			// maximum (stability constraint)
			if (irrev -> viscosity_curl_rhombi[vector_index] > irrev -> max_diff_h_coeff_turb)
			{
				irrev -> viscosity_curl_rhombi[vector_index] = irrev -> max_diff_h_coeff_turb;
			}
			
			// multiplying by the mass density of the gas phase
			irrev -> viscosity_curl_rhombi[vector_index] = 0.5*(density_gas(state, scalar_index_from) + density_gas(state, scalar_index_to))
			*irrev -> viscosity_curl_rhombi[vector_index] ;
		}
	}
	// averaging the curl diffusion coefficient from edges to cells
	edges_to_cells(irrev -> viscosity_curl_rhombi, irrev -> viscosity_curl, grid);
	return 0;
}

int hor_curl_viscosity_triangles(State *state, Irreversible_quantities *irrev, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Config *config)
{
	/*
	This function computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal curl movements on triangles.
	*/
	
	int layer_index, h_index, rho_base_index, scalar_base_index;
	double molecular_viscosity, density_value;
	#pragma omp parallel for private(molecular_viscosity, layer_index, h_index, density_value, rho_base_index, scalar_base_index)
	for (int i = 0; i < NO_OF_DUAL_V_VECTORS; ++i)
	{
		layer_index = i/NO_OF_DUAL_SCALARS_H;
		h_index = i - layer_index*NO_OF_DUAL_SCALARS_H;
		
		scalar_base_index = layer_index*NO_OF_SCALARS_H;
		
		// preliminary result
		irrev -> viscosity_curl_triangles[i] = 1.0/6.0*(
		tke2hor_diff_coeff(irrev -> tke[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]])
		+ tke2hor_diff_coeff(irrev -> tke[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]])
		+ tke2hor_diff_coeff(irrev -> tke[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]])
		+ tke2hor_diff_coeff(irrev -> tke[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]])
		+ tke2hor_diff_coeff(irrev -> tke[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]])
		+ tke2hor_diff_coeff(irrev -> tke[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]));
		
		// calculating and adding the molecular viscosity
		rho_base_index = NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H;
		density_value =
		1.0/6.0*(
		state -> rho[rho_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ state -> rho[rho_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ state -> rho[rho_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ state -> rho[rho_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ state -> rho[rho_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ state -> rho[rho_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]);
		
		molecular_viscosity = 1.0/6.0*(
		irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]);
		irrev -> viscosity_curl_triangles[i] += molecular_viscosity;
		
		// maximum (stability constraint)
		if (irrev -> viscosity_curl_triangles[i] > irrev -> max_diff_h_coeff_turb)
		{
			irrev -> viscosity_curl_triangles[i] = irrev -> max_diff_h_coeff_turb;
		}
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_curl_triangles[i] = density_value*irrev -> viscosity_curl_triangles[i];
	}
	return 0;
}

int vert_hor_mom_viscosity(State *state, Irreversible_quantities *irrev, Diagnostics *diagnostics, Config *config, Grid *grid, double delta_t)
{
	/*
	This function computes the effective viscosity (Eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	This quantity is located at the half level edges.
	To obey the symmetry of the stress tensor, the same coefficient must be used for the horizontal diffusion of vertical velocity.
	*/
	
	double max_diff_v_coeff_turb = 0.125*pow(
	grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H]
	, 2)/delta_t;
	int layer_index, h_index, scalar_base_index;
	double mom_diff_coeff, molecular_viscosity;
	// updating the TKE
	tke_update(irrev, delta_t, state, diagnostics, grid);
	// loop over horizontal vector points at half levels
	#pragma omp parallel for private(layer_index, h_index, mom_diff_coeff, molecular_viscosity, scalar_base_index)
	for (int i = 0; i < NO_OF_H_VECTORS - NO_OF_VECTORS_H; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		scalar_base_index = layer_index*NO_OF_SCALARS_H;
		// the turbulent component
		mom_diff_coeff = 0.25*(tke2vert_diff_coeff(irrev -> tke[scalar_base_index + grid -> from_index[h_index]])
		+ tke2vert_diff_coeff(irrev -> tke[scalar_base_index + grid -> to_index[h_index]])
		+ tke2vert_diff_coeff(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]])
		+ tke2vert_diff_coeff(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]));
		// computing and adding the molecular viscosity
		// the scalar variables need to be averaged to the vector points at half levels
		molecular_viscosity = 0.25*(irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> from_index[h_index]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> to_index[h_index]]
		+ irrev -> molecular_diffusion_coeff[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ irrev -> molecular_diffusion_coeff[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]);
		mom_diff_coeff += molecular_viscosity;
		
		// obeying the stability limit
		if (mom_diff_coeff > max_diff_v_coeff_turb)
		{
			mom_diff_coeff = max_diff_v_coeff_turb;
		}
		
		// multiplying by the density (averaged to the half level edge)
		irrev -> vert_hor_viscosity[i + NO_OF_VECTORS_H] = 
		0.25*(density_gas(state, scalar_base_index + grid -> from_index[h_index])
		+ density_gas(state, scalar_base_index + grid -> to_index[h_index])
		+ density_gas(state, (layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index])
		+ density_gas(state, (layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]))
		*mom_diff_coeff;
	}
	// for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		irrev -> vert_hor_viscosity[i] = irrev -> vert_hor_viscosity[i + NO_OF_VECTORS_H];
	}
	// for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
	#pragma omp parallel for	
	for (int i = NO_OF_H_VECTORS; i < NO_OF_H_VECTORS + NO_OF_VECTORS_H; ++i)
	{
		irrev -> vert_hor_viscosity[i] = irrev -> vert_hor_viscosity[i - NO_OF_VECTORS_H];
	}
	return 0;
}

int vert_vert_mom_viscosity(State *state, Grid *grid, Diagnostics *diagnostics, Irreversible_quantities *irrev, double delta_t)
{
	/*
	This function multiplies scalar_field_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
	*/
	// the maximum vertical diffusion coefficient
	double max_diff_v_coeff_turb = 0.125*pow(
	grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H]
	, 2)/delta_t;
	int i;
	double mom_diff_coeff;
	#pragma omp parallel for private(mom_diff_coeff, i)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			mom_diff_coeff
			// molecular viscosity
			= irrev -> molecular_diffusion_coeff[i]
			// turbulent component
			+ tke2vert_diff_coeff(irrev -> tke[i]);
			// stability criterion
			if (mom_diff_coeff > max_diff_v_coeff_turb)
			{
				mom_diff_coeff = max_diff_v_coeff_turb;
			}
			
			diagnostics -> scalar_field_placeholder[i] = density_gas(state, i)*mom_diff_coeff*diagnostics -> scalar_field_placeholder[i];
		}
	}
	return 0;
}

int temp_diffusion_coeffs(State *state, Config *config, Irreversible_quantities *irrev, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including eddies).
	*/
	// The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
	if (config -> momentum_diff_h == 0)
	{
		hor_div_viscosity(state, irrev, grid, diagnostics, config);
		hor_curl_viscosity_rhombi(state, irrev, grid, diagnostics, config);
		tke_update(irrev, delta_t, state, diagnostics, grid);
		// molecular viscosity
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irrev -> molecular_diffusion_coeff[i] = calc_diffusion_coeff(diagnostics -> temperature_gas[i],
			state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
		}
	}
	double c_g_v;
	#pragma omp parallel for private(c_g_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		c_g_v = spec_heat_cap_diagnostics_v(state, i, config);
		// horizontal diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_h[i]
		= c_g_v*0.5*(irrev -> viscosity_div[i] + irrev -> viscosity_curl[i]);
		// vertical diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_v[i]
		// molecular component
		= density_gas(state, i)*c_g_v*(irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ tke2vert_diff_coeff(irrev -> tke[i]));
	}
	return 0;
}

int mass_diffusion_coeffs(State *state, Config *config, Irreversible_quantities *irrev, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous tracer diffusion coefficient (including eddies).
	*/
	// The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
	if (config -> momentum_diff_h == 0 && config -> temperature_diff_h == 0)
	{
		hor_div_viscosity(state, irrev, grid, diagnostics, config);
		hor_curl_viscosity_rhombi(state, irrev, grid, diagnostics, config);
		tke_update(irrev, delta_t, state, diagnostics, grid);
		// molecular viscosity
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irrev -> molecular_diffusion_coeff[i] = calc_diffusion_coeff(diagnostics -> temperature_gas[i],
			state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// horizontal diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_h[i]
		= 0.5*(irrev -> viscosity_div[i] + irrev -> viscosity_curl[i])/density_gas(state, i);
		// vertical diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_v[i]
		// molecular component
		= irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ tke2vert_diff_coeff(irrev -> tke[i]);
	}
	return 0;
}

double tke2hor_diff_coeff(double tke)
{
	/*
	This function returns the horizontal kinematic Eddy viscosity as a function of the specific TKE.
	*/
	
	double prop_constant = 48400.0; // unit: m
	double result = prop_constant*pow(tke, 0.5);
	return result;
}

double tke2vert_diff_coeff(double tke)
{
	/*
	This function returns the vertical kinematic Eddy viscosity as a function of the specific TKE.
	*/
	
	double prop_constant = 0.4; // unit: m
	double result = prop_constant*pow(tke, 0.5);
	return result;
}






