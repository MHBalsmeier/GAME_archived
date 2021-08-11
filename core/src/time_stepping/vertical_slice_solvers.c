/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This file contains the implicit vertical solvers.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../thermodynamics/thermodynamics.h"
#include "../spatial_operators/spatial_operators.h"
#include "atmostracers.h"

int thomas_algorithm(double [], double [], double [], double [], double [], int);

int three_band_solver_ver_waves(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, Grid *grid, int rk_step)
{
	/*
	This is the implicit vertical solver for the main fluid constituent.
	*/
	
	// declaring and defining some variables that will be needed later on
	int upper_index, lower_index;
	double impl_weight = get_impl_thermo_weight();
	double c_v = spec_heat_capacities_v_gas(0);
	double c_p = spec_heat_capacities_p_gas(0);
	double r_d = specific_gas_constants(0);
	// This is for Klemp (2008).
	double damping_coeff, damping_start_height, z_above_damping;
	damping_start_height = config_info -> damping_start_height_over_toa*grid -> z_vector[0];
	
	// loop over all columns
	#pragma omp parallel for private(upper_index, lower_index, damping_coeff, z_above_damping)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the Kompendium
		double c_vector[NO_OF_LAYERS - 2];
		double d_vector[NO_OF_LAYERS - 1];
		double e_vector[NO_OF_LAYERS - 2];
		double r_vector[NO_OF_LAYERS - 1];
		double rho_expl[NO_OF_LAYERS];
		double rhotheta_expl[NO_OF_LAYERS];
		double theta_pert_expl[NO_OF_LAYERS];
		double exner_pert_expl[NO_OF_LAYERS];
		double theta_int_new[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS - 1];
		double rho_int_old[NO_OF_LAYERS - 1];
		double rho_int_expl[NO_OF_LAYERS - 1];
		double theta_int_expl[NO_OF_LAYERS - 1];
		double alpha_old[NO_OF_LAYERS];
		double beta_old[NO_OF_LAYERS];
		double gamma_old[NO_OF_LAYERS];
		double alpha_new[NO_OF_LAYERS];
		double beta_new[NO_OF_LAYERS];
		double gamma_new[NO_OF_LAYERS];
		double alpha[NO_OF_LAYERS];
		double beta[NO_OF_LAYERS];
		double gamma[NO_OF_LAYERS];
		double density_interface_new;
		
		// explicit quantities
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			// explicit density
			rho_expl[j] = state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			+ delta_t*state_tendency -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			// explicit potential temperature density
			rhotheta_expl[j] = state_old -> rhotheta[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> rhotheta[j*NO_OF_SCALARS_H + i];
			if (i == 0)
			{
				// old time step partial derivatives of rho*theta and Pi (divided by the volume)
				alpha[j] = -state_old -> rhotheta[i + j*NO_OF_SCALARS_H]/pow(state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H], 2)
				/grid -> volume[i + j*NO_OF_SCALARS_H];
				beta[j] = 1/state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H]/grid -> volume[i + j*NO_OF_SCALARS_H];
				gamma[j] = r_d/(c_v*state_old -> rhotheta[i + j*NO_OF_SCALARS_H])
				*(grid -> exner_bg[i + j*NO_OF_SCALARS_H] + state_old -> exner_pert[i + j*NO_OF_SCALARS_H])/grid -> volume[i + j*NO_OF_SCALARS_H];
			}
			else
			{
				// old time step partial derivatives of rho*theta and Pi
				alpha_old[j] = -state_old -> rhotheta[i + j*NO_OF_SCALARS_H]/pow(state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H], 2);
				beta_old[j] = 1/state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H];
				gamma_old[j] = r_d/(c_v*state_old -> rhotheta[i + j*NO_OF_SCALARS_H])*(grid -> exner_bg[i + j*NO_OF_SCALARS_H] + state_old -> exner_pert[i + j*NO_OF_SCALARS_H]);
				// new time step partial derivatives of rho*theta and Pi
				alpha_new[j] = -state_new -> rhotheta[i + j*NO_OF_SCALARS_H]/pow(state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H], 2);
				beta_new[j] = 1/state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H];
				gamma_new[j] = r_d/(c_v*state_new -> rhotheta[i + j*NO_OF_SCALARS_H])*(grid -> exner_bg[i + j*NO_OF_SCALARS_H] + state_new -> exner_pert[i + j*NO_OF_SCALARS_H]);
				// interpolation in time and dividing by the volume
				alpha[j] = ((1 - impl_weight)*alpha_old[j] + impl_weight*alpha_new[j])/grid -> volume[i + j*NO_OF_SCALARS_H];
				beta[j] = ((1 - impl_weight)*beta_old[j] + impl_weight*beta_new[j])/grid -> volume[i + j*NO_OF_SCALARS_H];
				gamma[j] = ((1 - impl_weight)*gamma_old[j] + impl_weight*gamma_new[j])/grid -> volume[i + j*NO_OF_SCALARS_H];
			}
			// explicit potential temperature perturbation
			theta_pert_expl[j] = state_old -> theta_pert[i + j*NO_OF_SCALARS_H] + delta_t*grid -> volume[i + j*NO_OF_SCALARS_H]*(
			alpha[j]*state_tendency -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + beta[j]*state_tendency -> rhotheta[j*NO_OF_SCALARS_H + i]);
			// explicit Exner pressure perturbation
			exner_pert_expl[j] = state_old -> exner_pert[i + j*NO_OF_SCALARS_H] + delta_t*grid -> volume[i + j*NO_OF_SCALARS_H]*gamma[j]*state_tendency -> rhotheta[j*NO_OF_SCALARS_H + i];
		}
		
		// determining the interface values
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			upper_index = i + j*NO_OF_SCALARS_H;
			lower_index = i + (j + 1)*NO_OF_SCALARS_H;
			rho_int_old[j] = 0.5*(state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + upper_index] + state_old -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + lower_index]);
			rho_int_expl[j] = 0.5*(rho_expl[j] + rho_expl[j + 1]);
			theta_int_expl[j] = 0.5*(grid -> theta_bg[upper_index] + theta_pert_expl[j] + grid -> theta_bg[lower_index] + theta_pert_expl[j + 1]);
			theta_int_new[j] = 0.5*(state_new -> rhotheta[upper_index]/state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + upper_index]
			+ state_new -> rhotheta[lower_index]/state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + lower_index]);
		}
		
		// filling up the coefficient vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{			
			// main diagonal
			d_vector[j] = -pow(theta_int_new[j], 2)*(gamma[j] + gamma[j + 1])
			+ 0.5*(grid -> exner_bg[j*NO_OF_SCALARS_H + i] - grid -> exner_bg[(j + 1)*NO_OF_SCALARS_H + i])
			*(alpha[j + 1] - alpha[j] + theta_int_new[j]*(beta[j + 1] - beta[j]))
			- (grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i])/(impl_weight*pow(delta_t, 2)*c_p*rho_int_old[j])
			*(2/grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]) - delta_t*state_old -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]*0.5
			*(1/grid -> volume[j*NO_OF_SCALARS_H + i] + 1/grid -> volume[(j + 1)*NO_OF_SCALARS_H + i]);
			// right hand side
			r_vector[j] = -(state_old -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i])
			*(grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i])
			/(impl_weight*pow(delta_t, 2)*c_p)
			+ theta_int_expl[j]*(exner_pert_expl[j] - exner_pert_expl[j + 1])/delta_t
			+ 0.5/delta_t*(theta_pert_expl[j] + theta_pert_expl[j + 1])*(grid -> exner_bg[j*NO_OF_SCALARS_H + i] - grid -> exner_bg[(j + 1)*NO_OF_SCALARS_H + i])
			- (grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i])/(impl_weight*pow(delta_t, 2)*c_p)
			*state_old -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]*rho_int_expl[j]/rho_int_old[j];
		}
		// Klemp (2008) upper boundary layer
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			z_above_damping = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - damping_start_height;
			if (z_above_damping < 0)
			{
				damping_coeff = 0;
			}
			else
			{
				damping_coeff = config_info -> damping_coeff_max*pow(sin(0.5*M_PI*z_above_damping/(grid -> z_vector[0] - damping_start_height)), 2);
			}
			d_vector[j] = d_vector[j]*(1 + delta_t*damping_coeff);
		}
		for (int j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			// lower diagonal
			c_vector[j] = theta_int_new[j + 1]*gamma[j + 1]*theta_int_new[j]
			+ 0.5*(grid -> exner_bg[(j + 1)*NO_OF_SCALARS_H + i] - grid -> exner_bg[(j + 2)*NO_OF_SCALARS_H + i])
			*(alpha[j + 1] + beta[j + 1]*theta_int_new[j])
			- (grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 2)*NO_OF_SCALARS_H + i])/(impl_weight*delta_t*c_p)*0.5
			*state_old -> wind[(j + 2)*NO_OF_VECTORS_PER_LAYER + i]/(grid -> volume[(j + 1)*NO_OF_SCALARS_H + i]*rho_int_old[j + 1]);
			// upper diagonal
			e_vector[j] = theta_int_new[j]*gamma[j + 1]*theta_int_new[j + 1]
			- 0.5*(grid -> exner_bg[j*NO_OF_SCALARS_H + i] - grid -> exner_bg[(j + 1)*NO_OF_SCALARS_H + i])
			*(alpha[j + 1] + beta[j + 1]*theta_int_new[j + 1])
			+ (grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i])/(impl_weight*delta_t*c_p)*0.5
			*state_old -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]/(grid -> volume[(j + 1)*NO_OF_SCALARS_H + i]*rho_int_old[j]);
		}
		
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_LAYERS - 1);
		
		/*
		Writing the result into the new state.
		--------------------------------------
		*/
		// mass density
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				= rho_expl[j] + delta_t*(solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				= rho_expl[j] + delta_t*(-solution_vector[j - 1])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else
			{
				state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				= rho_expl[j] + delta_t*(-solution_vector[j - 1] + solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
		}
		// potential temperature density
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				state_new -> rhotheta[j*NO_OF_SCALARS_H + i]
				= rhotheta_expl[j] + delta_t*(theta_int_new[j]*solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				state_new -> rhotheta[j*NO_OF_SCALARS_H + i]
				= rhotheta_expl[j] + delta_t*(-theta_int_new[j - 1]*solution_vector[j - 1])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else
			{
				state_new -> rhotheta[j*NO_OF_SCALARS_H + i]
				= rhotheta_expl[j] + delta_t*(-theta_int_new[j - 1]*solution_vector[j - 1] + theta_int_new[j]*solution_vector[j])
				/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
		}
		// vertical velocity
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			density_interface_new
			= 0.5*(state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			+ state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (j + 1)*NO_OF_SCALARS_H + i]);
			state_new -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]
			= (2*solution_vector[j]/grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - density_interface_new*state_old -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i])
			/rho_int_old[j];
		}
		// Exner pressure perturbation
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> exner_pert[j*NO_OF_SCALARS_H + i] = theta_pert_expl[j] + grid -> volume[j*NO_OF_SCALARS_H + i]
			*gamma[j]*(state_new -> rhotheta[j*NO_OF_SCALARS_H + i] - rhotheta_expl[j]);
		}
		// potential temperature perturbation
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> theta_pert[j*NO_OF_SCALARS_H + i] = state_new -> rhotheta[j*NO_OF_SCALARS_H + i]
			/state_new -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			- grid -> theta_bg[j*NO_OF_SCALARS_H + i];
		}
		
	} // end of the column (index i) loop
	return 0;
}

int three_band_solver_gen_densitites(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, Grid *grid)
{
	// Vertical advection of generalized densities (of tracers) with 3-band matrices.
	// mass densities, density x temperatures
	int no_of_relevant_constituents;
	double impl_weight, expl_weight;
	impl_weight = 0.5;
	expl_weight = 1 - impl_weight;
	for (int quantity_id = 0; quantity_id < 2; ++quantity_id)
	{
		no_of_relevant_constituents = 0;
		// mass densities
		if (quantity_id == 0)
		{
			// all constituents have a mass density
			no_of_relevant_constituents = NO_OF_CONSTITUENTS - 1;
		}
		// density x temperature fields
		if (quantity_id == 1)
		{
			// in this case, all the condensed constituents have a density x temperature field
			if(config_info -> assume_lte == 0)
			{
				no_of_relevant_constituents = NO_OF_CONDENSED_CONSTITUENTS;
			}
			// in this case, no density x temperature fields are taken into account
			else
			{
				no_of_relevant_constituents = 0;
			}
		}
		
		// loop over all relevant constituents
		for (int k = 0; k < no_of_relevant_constituents; ++k)
		{
			// This is done for all tracers apart from the main gaseous constituent.
		 	if (quantity_id != 0 || k != NO_OF_CONDENSED_CONSTITUENTS)
		 	{
				// loop over all columns
				#pragma omp parallel for
				for (int i = 0; i < NO_OF_SCALARS_H; ++i)
				{
					// for meanings of these vectors look into the definition of the function thomas_algorithm
					double c_vector[NO_OF_LAYERS - 1];
					double d_vector[NO_OF_LAYERS];
					double e_vector[NO_OF_LAYERS - 1];
					double r_vector[NO_OF_LAYERS];
					double vertical_flux_vector_impl[NO_OF_LAYERS - 1];
					double vertical_flux_vector_rhs[NO_OF_LAYERS - 1];
					double solution_vector[NO_OF_LAYERS];
					double density_old_at_interface, area;
					int lower_index, upper_index;
					
					// diagnozing the vertical fluxes
					for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
					{
						vertical_flux_vector_impl[j] = state_old -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
						vertical_flux_vector_rhs[j] = state_new -> wind[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
						// preparing the vertical interpolation
						lower_index = i + (j + 1)*NO_OF_SCALARS_H;
						upper_index = i + j*NO_OF_SCALARS_H;
						// For condensed constituents, a sink velocity must be added.
						if (k < NO_OF_CONDENSED_CONSTITUENTS)
						{
							vertical_flux_vector_impl[j] -= 1;
							vertical_flux_vector_rhs[j] -= 1;
						}
						// multiplying the vertical velocity by the area
						area = grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
						vertical_flux_vector_impl[j] = area*vertical_flux_vector_impl[j];
						vertical_flux_vector_rhs[j] = area*vertical_flux_vector_rhs[j];
						// old density at the interface
						density_old_at_interface
						= 0.5*(state_old -> rho[k*NO_OF_SCALARS + upper_index]
						+ state_old -> rho[k*NO_OF_SCALARS + lower_index]);
						vertical_flux_vector_rhs[j] = density_old_at_interface*vertical_flux_vector_rhs[j];
					}
					
					/*
					Now we proceed to solve the vertical tridiagonal problems.
					*/
					// filling up the original vectors
					for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
					{
						c_vector[j] = impl_weight*0.5*delta_t/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
						e_vector[j] = -impl_weight*0.5*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
					}
					for (int j = 0; j < NO_OF_LAYERS; ++j)
					{
						if (j == 0)
						{
							d_vector[j] = 1
							- impl_weight*0.5*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[0];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							d_vector[j] = 1
							+ impl_weight*0.5*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j - 1];
						}
						else
						{
							d_vector[j] = 1
							+ impl_weight*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]
							*0.5*(vertical_flux_vector_impl[j - 1] - vertical_flux_vector_impl[j]);
						}
						// the explicit component
						// mass densities
						if (quantity_id == 0)
						{
							r_vector[j] =
							state_old -> rho[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
							+ delta_t*state_tendency -> rho[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
						}
						// density x temperatures
						if (quantity_id == 1)
						{
							r_vector[j] =
							state_old -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
							+ delta_t*state_tendency -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
						}
						// adding the explicit part of the vertical flux divergence
						if (j == 0)
						{
							r_vector[j] += expl_weight*delta_t*vertical_flux_vector_rhs[j]/grid -> volume[j*NO_OF_SCALARS_H + i];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							r_vector[j] += -expl_weight*delta_t*vertical_flux_vector_rhs[j - 1]/grid -> volume[j*NO_OF_SCALARS_H + i];
						}
						else
						{
							r_vector[j] += expl_weight*delta_t*(-vertical_flux_vector_rhs[j - 1] + vertical_flux_vector_rhs[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
						}
					}
					
					// calling the algorithm to solve the system of linear equations
					thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_LAYERS);
					
					// writing the result into the new state
					for (int j = 0; j < NO_OF_LAYERS; ++j)
					{
						// limiter: none of the densities may be negative
						if (solution_vector[j] < 0)
						{
							solution_vector[j] = 0;
						}
						
						// mass densities
						if (quantity_id == 0)
						{
							state_new -> rho[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
						}
						
						// density x temperature fields
						if (quantity_id == 1)
						{
							state_new -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
						}
					}
				} // horizontal index
			}
		} // constituent
	} // quantity
	return 0;
}

int thomas_algorithm(double c_vector[], double d_vector[], double e_vector[], double r_vector[], double solution_vector[], int solution_length)
{
	/*
	This function solves a system of linear equations with a three-band matrix.
	*/
	
	double e_prime_vector[solution_length - 1];
	double r_prime_vector[solution_length];
	// downward sweep (matrix)
	if (d_vector[0] != 0)
	{
		e_prime_vector[0] = e_vector[0]/d_vector[0];
	}
	else
	{
		e_prime_vector[0] = 0;
	}
	for (int j = 1; j < solution_length - 1; ++j)
	{
		if (d_vector[j] - e_prime_vector[j - 1]*c_vector[j - 1] != 0)
		{
			e_prime_vector[j] = e_vector[j]/(d_vector[j] - e_prime_vector[j - 1]*c_vector[j - 1]);
		}
		else
		{
			e_prime_vector[j] = 0;
		}
	}
	// downward sweep (right-hand side)
	if (d_vector[0] != 0)
	{
		r_prime_vector[0] = r_vector[0]/d_vector[0];
	}
	else
	{
		r_prime_vector[0] = 0;
	}
	for (int j = 1; j < solution_length; ++j)
	{
		if (d_vector[j] - e_prime_vector[j - 1]*c_vector[j - 1] != 0)
		{
			r_prime_vector[j] = (r_vector[j] - r_prime_vector[j - 1]*c_vector[j - 1])/(d_vector[j] - e_prime_vector[j - 1]*c_vector[j - 1]);
		}
		else
		{
			r_prime_vector[j] = 0;
		}
	}
	
	// upward sweep (final solution)
	solution_vector[solution_length - 1] = r_prime_vector[solution_length - 1];
	for (int j = solution_length - 2; j >= 0; --j)
	{
		solution_vector[j] = r_prime_vector[j] - e_prime_vector[j]*solution_vector[j + 1];
	}
	return 0;
}



