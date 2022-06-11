/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, diffusion coefficients, including eddy viscosities, are computed.
*/

#include <math.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "subgrid_scale.h"

double tke2hor_diff_coeff(double, double);
double tke2vert_diff_coeff(double, double, double);

int hor_viscosity(State *state, Irreversible_quantities *irrev, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Config *config)
{
	/*
	This function computes the effective diffusion coefficient (molecular + turbulent).
	*/
			
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// molecular component
		irrev -> molecular_diffusion_coeff[i] = calc_diffusion_coeff(diagnostics -> temperature[i], state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
		irrev -> viscosity[i] = irrev -> molecular_diffusion_coeff[i];
		
		// computing and adding the turbulent component
		irrev -> viscosity[i] += tke2hor_diff_coeff(irrev -> tke[i], grid -> eff_hor_res);
	}
	
	/*
	Averaging the viscosity to rhombi
	---------------------------------
	*/
	int scalar_index_from, scalar_index_to, vector_index;
	#pragma omp parallel for private(scalar_index_from, scalar_index_to, vector_index)
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
			
			// indices of the adjacent scalar grid points
			scalar_index_from = layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index];
			scalar_index_to = layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index];
			
			// preliminary result
			irrev -> viscosity_rhombi[vector_index] = 0.5*(irrev -> viscosity[scalar_index_from] + irrev -> viscosity[scalar_index_to]);
			
			// multiplying by the mass density of the gas phase
			irrev -> viscosity_rhombi[vector_index] = 0.5*(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_from]
			+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index_to])
			*irrev -> viscosity_rhombi[vector_index] ;
		}
	}
	
	/*
	Averaging the viscosity to triangles
	------------------------------------
	*/
	int layer_index, h_index, rho_base_index, scalar_base_index;
	double density_value;
	#pragma omp parallel for private(layer_index, h_index, density_value, rho_base_index, scalar_base_index)
	for (int i = 0; i < NO_OF_DUAL_V_VECTORS; ++i)
	{
		layer_index = i/NO_OF_DUAL_SCALARS_H;
		h_index = i - layer_index*NO_OF_DUAL_SCALARS_H;
		
		scalar_base_index = layer_index*NO_OF_SCALARS_H;
		
		// preliminary result
		irrev -> viscosity_triangles[i] = 1.0/6.0*(
		irrev -> viscosity[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ irrev -> viscosity[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 0]]]
		+ irrev -> viscosity[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ irrev -> viscosity[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 1]]]
		+ irrev -> viscosity[scalar_base_index + grid -> from_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]
		+ irrev -> viscosity[scalar_base_index + grid -> to_index[dualgrid -> vorticity_indices_triangles[3*h_index + 2]]]);
		
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
		
		// multiplying by the mass density of the gas phase
		irrev -> viscosity_triangles[i] = density_value*irrev -> viscosity_triangles[i];
	}
	
	/*
	Multiplying the viscosity in the cell centers by the gas density
	----------------------------------------------------------------
	*/
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// multiplying by the density
		irrev -> viscosity[i] = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*tke2hor_diff_coeff(irrev -> tke[i], grid -> eff_hor_res);
	}
	
	return 0;
}

int vert_hor_mom_viscosity(State *state, Irreversible_quantities *irrev, Diagnostics *diagnostics, Config *config, Grid *grid, double delta_t)
{
	/*
	This function computes the effective viscosity (eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	This quantity is located at the half level edges.
	*/
	
	int layer_index, h_index, scalar_base_index;
	double mom_diff_coeff, molecular_viscosity;
	// loop over horizontal vector points at half levels
	#pragma omp parallel for private(layer_index, h_index, mom_diff_coeff, molecular_viscosity, scalar_base_index)
	for (int i = 0; i < NO_OF_H_VECTORS - NO_OF_VECTORS_H; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		scalar_base_index = layer_index*NO_OF_SCALARS_H;
		// the turbulent component
		mom_diff_coeff = 0.25*(tke2vert_diff_coeff(irrev -> tke[scalar_base_index + grid -> from_index[h_index]],
		diagnostics -> n_squared[scalar_base_index + grid -> from_index[h_index]], grid -> layer_thickness[scalar_base_index + grid -> from_index[h_index]])
		+ tke2vert_diff_coeff(irrev -> tke[scalar_base_index + grid -> to_index[h_index]],
		diagnostics -> n_squared[scalar_base_index + grid -> to_index[h_index]], grid -> layer_thickness[scalar_base_index + grid -> to_index[h_index]])
		+ tke2vert_diff_coeff(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]],
		diagnostics -> n_squared[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]], grid -> layer_thickness[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]])
		+ tke2vert_diff_coeff(irrev -> tke[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]],
		diagnostics -> n_squared[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]], grid -> layer_thickness[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]));
		// computing and adding the molecular viscosity
		// the scalar variables need to be averaged to the vector points at half levels
		molecular_viscosity = 0.25*(irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> from_index[h_index]]
		+ irrev -> molecular_diffusion_coeff[scalar_base_index + grid -> to_index[h_index]]
		+ irrev -> molecular_diffusion_coeff[(layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ irrev -> molecular_diffusion_coeff[(layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]);
		mom_diff_coeff += molecular_viscosity;
		
		// multiplying by the density (averaged to the half level edge)
		irrev -> vert_hor_viscosity[i + NO_OF_VECTORS_H] = 
		0.25*(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_base_index + grid -> from_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_base_index + grid -> to_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (layer_index + 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
		+ state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (layer_index + 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]])
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
			+ tke2vert_diff_coeff(irrev -> tke[i], diagnostics -> n_squared[i], grid -> layer_thickness[i]);
			
			diagnostics -> scalar_field_placeholder[i] = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*mom_diff_coeff*diagnostics -> scalar_field_placeholder[i];
		}
	}
	return 0;
}

int scalar_diffusion_coeffs(State *state, Config *config, Irreversible_quantities *irrev, Diagnostics *diagnostics, double delta_t, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This function computes the scalar diffusion coefficients (including eddies).
	*/
	
	// The diffusion coefficient only has to be calculated if it has not yet been done.
	if (config -> momentum_diff_h == 0)
	{
		hor_viscosity(state, irrev, grid, dualgrid, diagnostics, config);
	}
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		/*
		Computing the mass diffusion coefficient
		----------------------------------------
		*/
		// horizontal diffusion coefficient
		irrev -> mass_diffusion_coeff_numerical_h[i]
		= irrev -> viscosity[i]/state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
		// vertical diffusion coefficient
		irrev -> mass_diffusion_coeff_numerical_v[i]
		// molecular component
		= irrev -> molecular_diffusion_coeff[i]
		// turbulent component
		+ tke2vert_diff_coeff(irrev -> tke[i], diagnostics -> n_squared[i], grid -> layer_thickness[i]);
		
		/*
		Computing the temperature diffusion coefficient
		-----------------------------------------------
		*/
		irrev -> temp_diffusion_coeff_numerical_h[i] = c_v_mass_weighted_air(state, diagnostics, i)*irrev -> mass_diffusion_coeff_numerical_h[i];
		irrev -> temp_diffusion_coeff_numerical_v[i] = c_v_mass_weighted_air(state, diagnostics, i)*irrev -> mass_diffusion_coeff_numerical_v[i];
	}
	return 0;
}

int update_n_squared(State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function calculates the Brunt-Väisälä frequency.
	*/
	
	// calculating the full virtual potential temperature
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = grid -> theta_v_bg[i] + state -> theta_v_pert[i];
	}
	// vertical gradient of the full virtual potential temperature
	grad_vert_cov(diagnostics -> scalar_field_placeholder, diagnostics -> vector_field_placeholder, grid);
	// calculating the inverse full virtual potential temperature
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = 1.0/diagnostics -> scalar_field_placeholder[i];
	}
	scalar_times_vector_v(diagnostics -> scalar_field_placeholder, diagnostics -> vector_field_placeholder, diagnostics -> vector_field_placeholder, grid);
	
	// multiplying by the gravity acceleration
    int layer_index, h_index, vector_index;
	#pragma omp parallel for private(layer_index, h_index, vector_index)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_V_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        vector_index = h_index + layer_index*NO_OF_VECTORS_PER_LAYER;
        diagnostics -> vector_field_placeholder[vector_index]
        = grid -> gravity_m[vector_index]*diagnostics -> vector_field_placeholder[vector_index];
    }
    
    // averaging vertically to the scalar points
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        if (layer_index == 0)
        {
        	diagnostics -> n_squared[i] = diagnostics -> vector_field_placeholder[NO_OF_VECTORS_PER_LAYER + i];
        }
        else if (layer_index == NO_OF_LAYERS - 1)
        {
        	diagnostics -> n_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + i]
        	= diagnostics -> vector_field_placeholder[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + i];
        }
    	else
    	{
        	diagnostics -> n_squared[i]
        	= grid -> inner_product_weights[8*i + 6]*diagnostics -> vector_field_placeholder[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]
        	+ grid -> inner_product_weights[8*i + 7]*diagnostics -> vector_field_placeholder[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
    	}
    }
    
	return 0;
}

double tke2hor_diff_coeff(double tke, double effective_resolution)
{
	/*
	This function returns the horizontal kinematic eddy viscosity as a function of the specific TKE.
	*/
	
	double mean_velocity = pow(2.0*tke, 0.5);
	double mean_free_path = effective_resolution/4.0;
	double result = 1.0/6.0*mean_free_path*mean_velocity;
	return result;
}

double tke2vert_diff_coeff(double tke, double n_squared, double layer_thickness)
{
	/*
	This function returns the vertical kinematic eddy viscosity as a function of the specific TKE and the Brunt-Väisälä frequency.
	*/
	
	// vertical component of the turbulent kinetic energy
	double tke_vert = 1e-2*tke;
	
	double mean_velocity = pow(2.0*tke_vert, 0.5);
	// used Brunt-Väisälä frequency
	double n_used = pow(fmax(n_squared, 1e-4), 0.5);
	double mean_free_path = pow(2.0*tke_vert, 0.5)/n_used;
	mean_free_path = fmin(mean_free_path, layer_thickness);
	double result = 1.0/6.0*mean_free_path*mean_velocity;
	return result;
}






