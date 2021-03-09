/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "diagnostics.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int calc_diffusion_coeff(double temperature, double particle_mass, double denstiy, double particle_radius, double *result);
int calc_shear(State *, Diagnostics *, Grid *);

int calc_mass_diffusion_coeffs(State *state, Config_info *config_info, Scalar_field mass_diffusion_coeff_numerical_h, Scalar_field mass_diffusion_coeff_numerical_v)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including Eddys).
	*/
	if (config_info -> mass_diff_h == 1 || config_info -> mass_diff_v == 1)
	{
		double mean_particle_mass = mean_particle_masses_gas(0);
		double eff_particle_radius = 130e-12;
		double mass_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v;
		#pragma omp parallel for private (mass_diffusion_coeff, upturn_for_scale_h, upturn_for_scale_v)
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
		    calc_diffusion_coeff(state -> temperature_gas[i], mean_particle_mass, state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], eff_particle_radius, &mass_diffusion_coeff);
		    if (config_info -> mass_diff_h == 1)
		    {
		    	upturn_for_scale_h = 1;
	    	}
		    else
		    {
		    	upturn_for_scale_h = 0;
	    	}
		    if (config_info -> mass_diff_v == 1)
		    {
		    	upturn_for_scale_v = 1;
	    	}
		    else
		    {
		    	upturn_for_scale_v = 0;
	    	}
		    mass_diffusion_coeff_numerical_h[i] = upturn_for_scale_h*mass_diffusion_coeff;
		    mass_diffusion_coeff_numerical_v[i] = upturn_for_scale_v*mass_diffusion_coeff;
		}
	}
	return 0;
}

int calc_temp_diffusion_coeffs(State *state, Config_info *config_info, Irreversible_quantities *irreversible_quantities, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	/*
	This function computes the viscous temperature diffusion coefficient (including Eddys).
	*/
	double rho_g, c_g_v;
	// The Eddy viscosity coefficient only has to be calculated if it has not yet been done.
	if (config_info -> momentum_diff_h == 1)
	{
		hori_viscosity_eff(state, irreversible_quantities -> viscosity_eff, grid, diagnostics, config_info, delta_t);
	}
	edges_to_cells(irreversible_quantities -> viscosity_eff, irreversible_quantities -> scalar_diffusion_coeff_numerical_h, grid);
	#pragma omp parallel for private (rho_g, c_g_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		rho_g = density_gas(state, i);
		c_g_v = spec_heat_cap_diagnostics_v(state, i, config_info);
	    irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i] = rho_g*c_g_v* irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i];
	    // vertical Eddy viscosity is about two orders of magnitude smaller
	    irreversible_quantities -> scalar_diffusion_coeff_numerical_v[i] = 0.01*irreversible_quantities -> scalar_diffusion_coeff_numerical_h[i];
	}
	return 0;
}

int hori_viscosity_eff(State *state, Vector_field viscosity_eff, Grid *grid, Diagnostics *diagnostics, Config_info *config_info, double delta_t)
{
	// these things are hardly ever modified
	double eff_particle_radius = 130e-12;
	double mean_particle_mass = mean_particle_masses_gas(0);
	double molecular_viscosity;
	// the maximum diffusion coefficient (based on stability, including a 30 % safety margin)
	double max_diff_h_coeff_turb = (1 - 0.3)*0.25*grid -> mean_area_edge/delta_t;
	// the minimum "background" diffusion coefficient
	double min_diff_h_coeff_turb = grid -> mean_area_edge*config_info -> diff_h_smag_fac*config_info -> shear_bg;
	double viscosity_value;
	// calculatuing the shear
	calc_shear(state, diagnostics, grid);
	int layer_index, h_index;
	#pragma omp parallel for private(molecular_viscosity, layer_index, h_index)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		
		// preliminary result
		viscosity_value =  grid -> mean_area_edge*config_info -> diff_h_smag_fac*diagnostics -> shear[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
		
		/*
		Checking if the calculated value is "allowed".
		*/
		// calculating the molecular viscosity
		calc_diffusion_coeff(
		0.5*(state -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ state -> temperature_gas[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]]),
		mean_particle_mass,
		0.5*(state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]]),
		eff_particle_radius,
		&molecular_viscosity);
		// the molecular viscosity is the absolute minimum
		if (viscosity_value < molecular_viscosity)
		{
			viscosity_value = molecular_viscosity;
		}
		// turbulent minimum
		if (viscosity_value < min_diff_h_coeff_turb)
		{
			viscosity_value = min_diff_h_coeff_turb;
		}
		// stability
		if (viscosity_value > max_diff_h_coeff_turb)
		{
			viscosity_value = max_diff_h_coeff_turb;
		}
		
		// result
		viscosity_eff[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index] = viscosity_value;
	}
	return 0;
}

int ver_viscosity_eff(State *state, Vector_field viscosity_eff, Grid *grid)
{
	int layer_index, no_of_edges, k;
	double vertical_mom_viscosity_over_hor_viscosity = 0.01;
	for (int i = 0; i < NO_OF_LEVELS; ++i)
	{
		#pragma omp parallel for private(layer_index, no_of_edges, k)
		for (int j = 0; j < NO_OF_SCALARS_H; ++j)
		{
			viscosity_eff[i*NO_OF_VECTORS_PER_LAYER + j] = 0;
			no_of_edges = 6;
			if (j < NO_OF_PENTAGONS)
			{
				no_of_edges = 5;
			}
			layer_index = i;
			if (i == NO_OF_LEVELS - 1)
			{
				layer_index = i - 1;
			}
			for (int k = 0; k < NO_OF_EDGES; ++k)
			{
				viscosity_eff[i*NO_OF_VECTORS_PER_LAYER + j] += vertical_mom_viscosity_over_hor_viscosity
				*1.0/no_of_edges
				*viscosity_eff[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*j + k]];
			}
		}
	}
}

int calc_shear(State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function calculates the shear of the horizontal wind field at vector points.
	*/
	// calculating eastward and northward winds at edge
	calc_uv_at_edge(state -> velocity_gas, diagnostics -> u_at_edge, diagnostics -> v_at_edge, grid);
	// averaging the wind components to the cell centers
	edges_to_cells(diagnostics -> u_at_edge, diagnostics -> u_at_cell, grid);
	edges_to_cells(diagnostics -> v_at_edge, diagnostics -> v_at_cell, grid);
	// computing the horizontal gradient of the eastward wind
	grad_hor(diagnostics -> u_at_cell, diagnostics -> u_at_cell_grad, grid);
	// computing the horizontal gradient of the northward wind
	grad_hor(diagnostics -> v_at_cell, diagnostics -> v_at_cell_grad, grid);
	int layer_index, h_index;
	double comp_orth, comp_tang, dudx, dudy, dvdx, dvdy;
	#pragma omp parallel for private(layer_index, h_index, dudx, dudy, dvdx, dvdy, comp_orth, comp_tang)
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		// only the shear of the horizontal wind field is diagnozed
		if (h_index >= NO_OF_SCALARS_H)
		{
			// diagnozing u quantities
			comp_orth = diagnostics -> u_at_cell_grad[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
			tangential_wind(diagnostics -> u_at_cell_grad, layer_index, h_index, &comp_tang, grid);
			passive_turn(comp_orth, comp_tang, -grid -> direction[h_index], &dudx, &dudy);
			// diagnozing v quantities
			comp_orth = diagnostics -> v_at_cell_grad[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
			tangential_wind(diagnostics -> v_at_cell_grad, layer_index, h_index, &comp_tang, grid);
			passive_turn(comp_orth, comp_tang, -grid -> direction[h_index], &dvdx, &dvdy);
			// calculating the deformation according to the MPAS paper
			diagnostics -> shear[i] = sqrt(pow(dudx - dvdy, 2) + pow(dudy + dvdx, 2));
		}
	}
	return 0;
}


int calc_diffusion_coeff(double temperature, double particle_mass, double density, double particle_radius, double *result)
{
    double thermal_velocity = sqrt(8*K_B*temperature/(M_PI*particle_mass));
    double particle_density = density/particle_mass;
    double cross_section = 4*M_PI*pow(particle_radius, 2);
    double mean_free_path = 1/(sqrt(2)*particle_density*cross_section);
    *result = 1.0/3.0*thermal_velocity*mean_free_path;
    return 0;
}







