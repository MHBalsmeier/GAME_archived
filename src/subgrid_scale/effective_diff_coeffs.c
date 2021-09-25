/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, diffusion coefficients, including Eddy viscosities, are computed.
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int tke_update(Irreversible_quantities *, double, State *, Diagnostics *, Grid *);
double return_ver_hor_viscosity(double);

int hori_div_viscosity_eff(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_cell*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	// the maximum diffusion coefficient (stability constraint)
	double max_diff_h_coeff_turb = 0.125*grid -> mean_area_cell/delta_t;
	double molecular_viscosity;
	
	#pragma omp parallel for private(molecular_viscosity)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// preliminary result
		irrev -> viscosity_div_eff[i] = 7*grid -> mean_area_cell*config_info -> diff_h_smag_fac
		*fabs(5.0/3*diagnostics -> wind_divv[i]);
		
		// calculating and adding the molecular viscosity
		calc_diffusion_coeff(diagnostics -> temperature_gas[i], mean_particle_mass, state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i],
		eff_particle_radius, &molecular_viscosity);
		irrev -> viscosity_div_eff[i] += molecular_viscosity;
		
		// turbulent minimum
		if (irrev -> viscosity_div_eff[i] < min_diff_h_coeff_turb)
		{
			irrev -> viscosity_div_eff[i] = min_diff_h_coeff_turb;
		}
		
		// maximum (stability constraint)
		if (irrev -> viscosity_div_eff[i] > max_diff_h_coeff_turb)
		{
			irrev -> viscosity_div_eff[i] = max_diff_h_coeff_turb;
		}
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_div_eff[i] = density_gas(state, i)*irrev -> viscosity_div_eff[i];
	}
	return 0;
}

int hori_curl_viscosity_eff_rhombi(State *state, Irreversible_quantities *irrev, Grid *grid, Diagnostics *diagnostics, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_cell*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	// the maximum diffusion coefficient (stability constraint)
	double max_diff_h_coeff_turb = 0.125*grid -> mean_area_cell/delta_t;
	double molecular_viscosity;
	
	int layer_index, h_index, scalar_index_from, scalar_index_to, vector_index;
	#pragma omp parallel for private(molecular_viscosity, layer_index, h_index, scalar_index_from, scalar_index_to, vector_index)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
		// preliminary result
		irrev -> viscosity_curl_eff_rhombi[vector_index] = 0.35*grid -> mean_area_cell*config_info -> diff_h_smag_fac
		*fabs(diagnostics -> rel_vort[NO_OF_VECTORS_H + 2*layer_index*NO_OF_VECTORS_H + h_index]);
		
		// calculating and adding the molecular viscosity
		scalar_index_from = layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index];
		scalar_index_to = layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index];
		calc_diffusion_coeff(
		0.5*(diagnostics -> temperature_gas[scalar_index_from] + diagnostics -> temperature_gas[scalar_index_to]),
		mean_particle_mass,
		0.5*(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_from] + state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_to]),
		eff_particle_radius, &molecular_viscosity);
		irrev -> viscosity_curl_eff_rhombi[vector_index] += molecular_viscosity;
		
		// turbulent minimum
		if (irrev -> viscosity_curl_eff_rhombi[vector_index] < min_diff_h_coeff_turb)
		{
			irrev -> viscosity_curl_eff_rhombi[vector_index] = min_diff_h_coeff_turb;
		}
		
		// maximum (stability constraint)
		if (irrev -> viscosity_curl_eff_rhombi[vector_index] > max_diff_h_coeff_turb)
		{
			irrev -> viscosity_curl_eff_rhombi[vector_index] = max_diff_h_coeff_turb;
		}
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_curl_eff_rhombi[vector_index] = 0.5*(density_gas(state, scalar_index_from) + density_gas(state, scalar_index_to))*irrev -> viscosity_curl_eff_rhombi[vector_index];
	}
	return 0;
}

int hori_curl_viscosity_eff_triangles(State *state, Irreversible_quantities *irrev, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_cell*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	// the maximum diffusion coefficient (stability constraint)
	double max_diff_h_coeff_turb = 0.125*grid -> mean_area_cell/delta_t;
	
	int layer_index, h_index;
	double molecular_viscosity, density_value;
	#pragma omp parallel for private(molecular_viscosity, layer_index, h_index, density_value)
	for (int i = 0; i < NO_OF_DUAL_V_VECTORS; ++i)
	{
		layer_index = i/NO_OF_DUAL_SCALARS_H;
		h_index = i - layer_index*NO_OF_DUAL_SCALARS_H;
		// preliminary result
		irrev -> viscosity_curl_eff_triangles[i] = 0.35*grid -> mean_area_cell*config_info -> diff_h_smag_fac
		*fabs(diagnostics -> rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + h_index]);
		
		// calculating and adding the molecular viscosity
		density_value =
		1.0/6*(
		state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]);
		calc_diffusion_coeff(
		1.0/6*(
		diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]),
		mean_particle_mass,
		density_value,
		eff_particle_radius, &molecular_viscosity);
		irrev -> viscosity_curl_eff_triangles[i] += molecular_viscosity;
		
		// turbulent minimum
		if (irrev -> viscosity_curl_eff_triangles[i] < min_diff_h_coeff_turb)
		{
			irrev -> viscosity_curl_eff_triangles[i] = min_diff_h_coeff_turb;
		}
		
		// maximum (stability constraint)
		if (irrev -> viscosity_curl_eff_triangles[i] > max_diff_h_coeff_turb)
		{
			irrev -> viscosity_curl_eff_triangles[i] = max_diff_h_coeff_turb;
		}
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_curl_eff_triangles[i] = density_value*irrev -> viscosity_curl_eff_triangles[i];
	}
	return 0;
}

int vert_hor_mom_viscosity(State *state, Irreversible_quantities *irrev, Diagnostics *diagnostics, Config_info *config_info, Grid *grid, double delta_t)
{
	/*
	This function computes the effective viscosity (Eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	This quantity is located at the half level edges.
	To obey the symmetry of the stress tensor, the same coefficient must be used for the horizontal diffusion of vertical velocity.
	*/
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	double max_diff_v_coeff_turb = 0.125*pow(
	grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H]
	, 2)/delta_t;
	int layer_index, h_index;
	double mom_diff_coeff, molecuar_viscosity;
	// updating the TKE
	tke_update(irrev, delta_t, state, diagnostics, grid);
	// loop over horizontal vector points at half levels
	#pragma omp parallel for private(layer_index, h_index, mom_diff_coeff, molecuar_viscosity)
	for (int i = 0; i < NO_OF_H_VECTORS - NO_OF_VECTORS_H; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		// the turbulent component
		mom_diff_coeff = 0.25*(return_ver_hor_viscosity(irrev -> tke[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]])
		+ return_ver_hor_viscosity(irrev -> tke[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]])
		+ return_ver_hor_viscosity(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]])
		+ return_ver_hor_viscosity(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]));
		// computing and adding the molecular viscosity
		// the scalar variables need to be averaged to the vector points at half levels
		calc_diffusion_coeff(0.25*(diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ diagnostics -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]]
		+ diagnostics -> temperature_gas[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ diagnostics -> temperature_gas[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]),
		mean_particle_mass,
		0.25*(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]), eff_particle_radius, &molecuar_viscosity);
		mom_diff_coeff += molecuar_viscosity;
		
		// obeying the stability limit
		if (mom_diff_coeff > max_diff_v_coeff_turb)
		{
			mom_diff_coeff = max_diff_v_coeff_turb;
		}
		
		// multiplying by the density (averaged to the half level edge)
		irrev -> vert_hor_viscosity_eff[i] = 
		0.25*(density_gas(state, layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index])
		+ density_gas(state, layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index])
		+ density_gas(state, (layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index])
		+ density_gas(state, (layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]))
		*mom_diff_coeff;
	}
	return 0;
}

int vert_w_viscosity_eff(State *state, Grid *grid, Diagnostics *diagnostics, double delta_t)
{
	/*
	This function multiplies scalar_field_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
	*/
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	// the maximum vertical diffusion coefficient
	double max_diff_v_coeff_turb = 0.125*pow(
	grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H]
	, 2)/delta_t;
	int layer_index, h_index;
	double mom_diff_coeff, molecuar_viscosity;
	#pragma omp parallel for private(mom_diff_coeff, molecuar_viscosity, h_index, layer_index)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		// this is the value resulting from turbulence
		mom_diff_coeff = 0.11*pow(
		grid -> z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]
		, 2)
		*fabs(diagnostics -> scalar_field_placeholder[i]);
		// computing and adding the molecular momentum diffusion
		calc_diffusion_coeff(diagnostics -> temperature_gas[i], mean_particle_mass,
		state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &molecuar_viscosity);
		mom_diff_coeff += molecuar_viscosity;
		// this is the stability criterion
		if (mom_diff_coeff > max_diff_v_coeff_turb)
		{
			mom_diff_coeff = max_diff_v_coeff_turb;
		}
		
		diagnostics -> scalar_field_placeholder[i] = density_gas(state, i)*mom_diff_coeff*diagnostics -> scalar_field_placeholder[i];
	}
	return 0;
}

int calc_temp_diffusion_coeffs(State *state, Config_info *config_info, Irreversible_quantities *irreversible_quantities, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including Eddys).
	*/
	double mean_particle_mass = mean_particle_masses_gas(0);
	double eff_particle_radius = 130e-12;
	// The Eddy viscosity coefficient only has to be calculated if it has not yet been done.
	if (config_info -> momentum_diff_h == 0)
	{
		hori_div_viscosity_eff(state, irreversible_quantities, grid, diagnostics, config_info, delta_t);
		hori_curl_viscosity_eff_rhombi(state, irreversible_quantities, grid, diagnostics, config_info, delta_t);
	}
	// averaging the curl diffusion coefficient from edges to cells
	edges_to_cells(irreversible_quantities -> viscosity_curl_eff_rhombi, diagnostics -> scalar_field_placeholder, grid);
	double c_g_v;
	#pragma omp parallel for private (c_g_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		c_g_v = spec_heat_cap_diagnostics_v(state, i, config_info);
		irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i] = c_g_v*(irreversible_quantities -> viscosity_div_eff[i] + diagnostics -> scalar_field_placeholder[i]);
		// the vertical viscosity is just the molecular viscosity for now
		calc_diffusion_coeff(diagnostics -> temperature_gas[i], mean_particle_mass,
		state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &irreversible_quantities -> scalar_diffusion_coeff_numerical_v[i]);
		irreversible_quantities -> scalar_diffusion_coeff_numerical_v[i] = c_g_v*irreversible_quantities -> scalar_diffusion_coeff_numerical_v[i];
	}
	return 0;
}

int calc_mass_diffusion_coeffs(State *state, Config_info *config_info, Irreversible_quantities *irreversible_quantities, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous tracer diffusion coefficient (including Eddys).
	*/
	// The Eddy viscosity coefficient only has to be calculated if it has not yet been done.
	if (config_info -> momentum_diff_h == 0)
	{
		hori_div_viscosity_eff(state, irreversible_quantities, grid, diagnostics, config_info, delta_t);
		hori_curl_viscosity_eff_rhombi(state, irreversible_quantities, grid, diagnostics, config_info, delta_t);
	}
	// averaging the curl diffusion coefficient from edges to cells
	edges_to_cells(irreversible_quantities -> viscosity_curl_eff_rhombi, diagnostics -> scalar_field_placeholder, grid);
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i] = 0.1*(irreversible_quantities -> viscosity_div_eff[i] + diagnostics -> scalar_field_placeholder[i])
		/density_gas(state, i);
		// the vertical viscosity is proportional to the horizontal viscosity for now
		irreversible_quantities -> scalar_diffusion_coeff_numerical_v[i] = 0.001*irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i];
	}
	return 0;
}

int tke_update(Irreversible_quantities *irrev, double delta_t, State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function updates the specific turbulent kinetic energy (TKE), unit: J/kg.
	*/
	// e-folding time of TKE over flat ground for RES_ID = 5
	double e_folding_time_flat = 86400.0;
	// e-folding time of TKE over rough ground for RES_ID = 5
	double e_folding_time_rough = 43200.0;
	// the decay time gets shorter for smaller mesh sizes
	double decay_constant_sea = 1.0/e_folding_time_flat*pow(2, RES_ID - 5);
	double decay_constant_land = 1.0/e_folding_time_rough*pow(2, RES_ID - 5);
	int i;
	double decay_constant;
	// computing the advection
	grad(irrev -> tke, diagnostics -> vector_field_placeholder, grid);
	inner_product(diagnostics -> vector_field_placeholder, state -> wind, diagnostics -> scalar_field_placeholder, grid);
	double production_rate;
	#pragma omp parallel for private(i, decay_constant, production_rate)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			production_rate = 0;
			// the decay constants differ over land vs over water
			if (grid -> is_land[h_index] == 1 && grid -> z_scalar[i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index] <= 1000.0)
			{
				decay_constant = decay_constant_land;
				production_rate = 0.5*decay_constant;
			}
			else
			{
				decay_constant = decay_constant_sea;
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
			+ production_rate*diagnostics -> e_kin[i]);
			if (irrev -> tke[i] < 0)
			{
				irrev -> tke[i] = 0;
			}
		}
	}
	return 0;
}

double return_ver_hor_viscosity(double tke)
{
	/*
	This function returns the vertical kinematic Eddy viscosity as a function of the specific TKE.
	*/
	double prop_constant = 32*pow(2, -RES_ID)*0.5; // unit: m
	double result = prop_constant*pow(tke, 0.5);
	return result;
}












