/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, diffusion coefficients, including Eddy viscosities, are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"

int tke_update(Irreversible_quantities *, double, State *, Diagnostics *, Grid *);
double ver_hor_viscosity(double, double, double);

int hori_div_viscosity(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config *config)
{
	/*
	This function computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal divergent movements.
	*/
	
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// molecular component
		irrev -> molecular_diffusion_coeff[i] = calc_diffusion_coeff(diagnostics -> temperature_gas[i], state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
		irrev -> viscosity_div[i] = irrev -> molecular_diffusion_coeff[i];
		
		// turbulent component (the divergence is approximately one order of magnitude smaller than the vorticity, that is where the prefactor of 6 comes from)
		// 4.0/3 is a result of the stress tensor
		irrev -> viscosity_div[i] += config -> diff_h_smag_div*grid -> mean_velocity_area*fabs(4.0/3*diagnostics -> wind_divv[i]);
		
		// maximum (stability constraint)
		if (irrev -> viscosity_div[i] > irrev -> max_diff_h_coeff_turb)
		{
			irrev -> viscosity_div[i] = irrev -> max_diff_h_coeff_turb;
		}
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_div[i] = density_gas(state, i)*irrev -> viscosity_div[i];
	}
	return 0;
}

int hori_curl_viscosity_rhombi(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config *config)
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
			// preliminary result
			irrev -> viscosity_curl_rhombi[vector_index] = config -> diff_h_smag_rot*grid -> mean_velocity_area
			*fabs(diagnostics -> rel_vort[NO_OF_VECTORS_H + 2*layer_index*NO_OF_VECTORS_H + h_index]);
			
			// calculating and adding the molecular viscosity
			scalar_index_from = layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index];
			scalar_index_to = layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index];
			molecular_viscosity = calc_diffusion_coeff(
			0.5*(diagnostics -> temperature_gas[scalar_index_from] + diagnostics -> temperature_gas[scalar_index_to]),
			0.5*(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_from] + state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_to]));
			irrev -> viscosity_curl_rhombi[vector_index] += molecular_viscosity;
			
			// maximum (stability constraint)
			if (irrev -> viscosity_curl_rhombi[vector_index] > irrev -> max_diff_h_coeff_turb)
			{
				irrev -> viscosity_curl_rhombi[vector_index] = irrev -> max_diff_h_coeff_turb;
			}
			
			// multiplying by the mass density of the gas phase
			irrev -> viscosity_curl_rhombi[vector_index] = 0.5*(density_gas(state, scalar_index_from) + density_gas(state, scalar_index_to))*irrev -> viscosity_curl_rhombi[vector_index];
		}
	}
	// averaging the curl diffusion coefficient from edges to cells
	edges_to_cells(irrev -> viscosity_curl_rhombi, irrev -> viscosity_curl, grid);
	return 0;
}

int hori_curl_viscosity_triangles(State *state, Irreversible_quantities *irrev, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Config *config)
{
	/*
	This function computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal curl movements on triangles.
	*/
	
	int layer_index, h_index, rho_base_index, temp_base_index;
	double molecular_viscosity, density_value;
	#pragma omp parallel for private(molecular_viscosity, layer_index, h_index, density_value, rho_base_index, temp_base_index)
	for (int i = 0; i < NO_OF_DUAL_V_VECTORS; ++i)
	{
		layer_index = i/NO_OF_DUAL_SCALARS_H;
		h_index = i - layer_index*NO_OF_DUAL_SCALARS_H;
		// preliminary result
		irrev -> viscosity_curl_triangles[i] = config -> diff_h_smag_rot*grid -> mean_velocity_area
		*fabs(diagnostics -> rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + h_index]);
		
		rho_base_index = NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H;
		// calculating and adding the molecular viscosity
		density_value =
		1.0/6*(
		state -> rho[rho_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ state -> rho[rho_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ state -> rho[rho_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ state -> rho[rho_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ state -> rho[rho_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ state -> rho[rho_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]);
		temp_base_index = layer_index*NO_OF_SCALARS_H ;
		molecular_viscosity = calc_diffusion_coeff(
		1.0/6*(
		diagnostics -> temperature_gas[temp_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ diagnostics -> temperature_gas[temp_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ diagnostics -> temperature_gas[temp_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ diagnostics -> temperature_gas[temp_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ diagnostics -> temperature_gas[temp_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ diagnostics -> temperature_gas[temp_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]),
		density_value);
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
	double mom_diff_coeff, molecuar_viscosity, delta_z;
	// updating the TKE
	tke_update(irrev, delta_t, state, diagnostics, grid);
	// loop over horizontal vector points at half levels
	#pragma omp parallel for private(layer_index, h_index, mom_diff_coeff, molecuar_viscosity, scalar_base_index, delta_z)
	for (int i = 0; i < NO_OF_H_VECTORS - NO_OF_VECTORS_H; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		scalar_base_index = layer_index*NO_OF_SCALARS_H;
		// the turbulent component
		delta_z = grid -> z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index]
		- grid -> z_vector[NO_OF_SCALARS_H + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER + h_index];
		mom_diff_coeff = 0.25*(ver_hor_viscosity(irrev -> tke[scalar_base_index + grid -> from_index[h_index]], delta_z, irrev -> mixing_length)
		+ ver_hor_viscosity(irrev -> tke[scalar_base_index + grid -> to_index[h_index]], delta_z, irrev -> mixing_length)
		+ ver_hor_viscosity(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]], delta_z, irrev -> mixing_length)
		+ ver_hor_viscosity(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]], delta_z, irrev -> mixing_length));
		// computing and adding the molecular viscosity
		// the scalar variables need to be averaged to the vector points at half levels
		molecuar_viscosity = calc_diffusion_coeff(0.25*(diagnostics -> temperature_gas[scalar_base_index + grid -> from_index[h_index]]
		+ diagnostics -> temperature_gas[scalar_base_index + grid -> to_index[h_index]]
		+ diagnostics -> temperature_gas[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ diagnostics -> temperature_gas[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]),
		0.25*(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_base_index + grid -> from_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_base_index + grid -> to_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]));
		mom_diff_coeff += molecuar_viscosity;
		
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

int vert_w_viscosity(State *state, Grid *grid, Diagnostics *diagnostics, Irreversible_quantities *irrev, double delta_t)
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
			+ pow(
			fmin(irrev -> mixing_length,
			grid -> z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER])
			, 2)
			*fabs(diagnostics -> scalar_field_placeholder[i]);
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

int calc_temp_diffusion_coeffs(State *state, Config *config, Irreversible_quantities *irrev, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including eddies).
	*/
	// The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
	if (config -> momentum_diff_h == 0)
	{
		hori_div_viscosity(state, irrev, grid, diagnostics, config);
		hori_curl_viscosity_rhombi(state, irrev, grid, diagnostics, config);
		tke_update(irrev, delta_t, state, diagnostics, grid);
		// molecular viscosity
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irrev -> molecular_diffusion_coeff[i] = calc_diffusion_coeff(diagnostics -> temperature_gas[i],
			state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
		}
	}
	int layer_index, h_index;
	double c_g_v;
	#pragma omp parallel for private(layer_index, h_index, c_g_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		c_g_v = spec_heat_cap_diagnostics_v(state, i, config);
		// horizontal diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_h[i]
		// molecular component
		= c_g_v*(density_gas(state, i)*irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ irrev -> viscosity_div[i] + irrev -> viscosity_curl[i]);
		// vertical diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_v[i]
		// molecular component
		= density_gas(state, i)*c_g_v*(irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ ver_hor_viscosity(irrev -> tke[i],
		grid -> z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER],
		irrev -> mixing_length));
	}
	return 0;
}

int calc_mass_diffusion_coeffs(State *state, Config *config, Irreversible_quantities *irrev, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous tracer diffusion coefficient (including eddies).
	*/
	// The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
	if (config -> momentum_diff_h == 0 && config -> temperature_diff_h == 0)
	{
		hori_div_viscosity(state, irrev, grid, diagnostics, config);
		hori_curl_viscosity_rhombi(state, irrev, grid, diagnostics, config);
		tke_update(irrev, delta_t, state, diagnostics, grid);
		// molecular viscosity
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			irrev -> molecular_diffusion_coeff[i] = calc_diffusion_coeff(diagnostics -> temperature_gas[i],
			state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
		}
	}
	int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		// horizontal diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_h[i]
		// molecular component
		= irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ (irrev -> viscosity_div[i] + irrev -> viscosity_curl[i])
		/density_gas(state, i);
		// vertical diffusion coefficient
		irrev -> scalar_diffusion_coeff_numerical_v[i]
		// molecular component
		= irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ ver_hor_viscosity(irrev -> tke[i], grid -> z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER],
		irrev -> mixing_length);
	}
	return 0;
}

int tke_update(Irreversible_quantities *irrev, double delta_t, State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function updates the specific turbulent kinetic energy (TKE), unit: J/kg.
	*/
	// the ratio of global unresolved to resolved kinetic energy
	double tke_ke_ratio = 0.1*pow(4, 5 - RES_ID);
	// the e-folding time of TKE approximation
	double tke_approx_time = 10800*pow(4, 5 - RES_ID);
	// computing the advection
	grad(irrev -> tke, diagnostics -> vector_field_placeholder, grid);
	inner_product(diagnostics -> vector_field_placeholder, state -> wind, diagnostics -> scalar_field_placeholder, grid);
	double boundary_layer_height = 1000.0;
	double roughness_length_factor = 1.0/0.08;
	int i;
	double decay_constant, production_rate;
	#pragma omp parallel for private(i, decay_constant, production_rate)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			
			// decay constant, as derived from diffusion
			decay_constant = 8*pow(M_PI, 2)/grid -> mean_velocity_area*(irrev -> viscosity_div[i] + irrev -> viscosity_curl[i])/density_gas(state, i);
			
			production_rate = 0;
			if (grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index] <= boundary_layer_height)
			{
				production_rate =
				// factor taking into account the roughness of the surface
				roughness_length_factor*grid -> roughness_length[h_index]
				// height-dependent factor
				*(boundary_layer_height - (grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index]))/boundary_layer_height
				*(tke_ke_ratio*0.5*diagnostics -> v_squared[i] - irrev -> tke[i])/tke_approx_time;
				// restricting the production rate to positive values
				production_rate = fmax(0, production_rate);
			}
			
			// prognostic equation for TKE
			irrev -> tke[i] += delta_t*(
			// advection
			- diagnostics -> scalar_field_placeholder[i]
			// production through dissipation of resolved energy
			+ irrev -> heating_diss[i]/density_gas(state, i)
			// decay through molecular dissipation
			- decay_constant*irrev -> tke[i]
			// production through turbulence generation in the boundary layer
			+ production_rate);
			// clipping negative values which might occur through advection
			if (irrev -> tke[i] < 0)
			{
				irrev -> tke[i] = 0;
			}
		}
	}
	return 0;
}

double ver_hor_viscosity(double tke, double delta_z, double mixing_length)
{
	/*
	This function returns the vertical kinematic Eddy viscosity as a function of the specific TKE.
	*/
	double prop_constant = 0.03*fmin(delta_z, mixing_length); // unit: m
	double result = prop_constant*pow(2*tke, 0.5);
	return result;
}












