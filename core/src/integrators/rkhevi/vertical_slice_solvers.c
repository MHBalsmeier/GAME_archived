/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This file contains the implicit vertical solvers.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../../settings.h"
#include "../../diagnostics/diagnostics.h"
#include "../../spatial_operators/spatial_operators.h"
#include "atmostracers.h"

int thomas_algorithm(double [], double [], double [], double [], double [], int);

int three_band_solver_ver_waves(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, Grid *grid)
{
	double upper_volume, lower_volume, total_volume, damping_coeff, damping_start_height, z_above_damping;
	// This is for Klemp (2008).
	damping_start_height = config_info -> damping_start_height_over_toa*grid -> z_vector[0];
	int upper_index, lower_index, j;
	double impl_pgrad_weight = get_impl_thermo_weight();
	double c_g_v = spec_heat_cap_diagnostics_v(state_old, NO_OF_SCALARS/2, config_info);
	double c_g_p = spec_heat_cap_diagnostics_p(state_old, NO_OF_SCALARS/2, config_info);
	#pragma omp parallel for private(upper_index, lower_index, j, upper_volume, lower_volume, total_volume, damping_coeff)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the Kompendium
		double c_vector[NO_OF_LAYERS - 2];
		double d_vector[NO_OF_LAYERS - 1];
		double e_vector[NO_OF_LAYERS - 2];
		double r_vector[NO_OF_LAYERS - 1];
		double density_explicit[NO_OF_LAYERS];
		double entropy_density_explicit[NO_OF_LAYERS];
		double spec_entropy_explicit[NO_OF_LAYERS];
		double density_new[NO_OF_LAYERS];
		double entropy_density_new[NO_OF_LAYERS];
		double spec_entropy_new[NO_OF_LAYERS];
		double spec_entropy_interface_new[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS - 1];
		double a[NO_OF_LAYERS - 1];
		double b[NO_OF_LAYERS - 1];
		double temp_old_interface_values[NO_OF_LAYERS - 1];
		double density_interface_old[NO_OF_LAYERS - 1];
		double density_interface_explicit[NO_OF_LAYERS - 1];
		double r_g_vector[NO_OF_LAYERS];
		double alpha[NO_OF_LAYERS];
		double beta[NO_OF_LAYERS];
		double gamma[NO_OF_LAYERS];
		double delta[NO_OF_LAYERS];
		double density_interface_new;
		// determining the properties of the gas phase in the grid boxes and explicit quantities
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			r_g_vector[j] = gas_constant_diagnostics(state_old, i + j*NO_OF_SCALARS_H, config_info);
			// explicit density
			density_explicit[j] = state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			+ delta_t*state_tendency -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			// explicit entropy density
			entropy_density_explicit[j] = state_old -> entropy_densities[j*NO_OF_SCALARS_H + i]
			+ delta_t*state_tendency -> entropy_densities[j*NO_OF_SCALARS_H + i];
			// specific entropy
			spec_entropy_explicit[j] = entropy_density_explicit[j]/density_explicit[j];
			// new density
			density_explicit[j] = state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			// new entropy density
			entropy_density_explicit[j] = state_new -> entropy_densities[j*NO_OF_SCALARS_H + i];
			// specific entropy
			spec_entropy_new[j] = entropy_density_new[j]/density_new[j];
			// partial derivatives of T = T(rho, stilde)
			alpha[j] = state_old -> temperature_gas[j*NO_OF_SCALARS_H + i]/(c_g_v*state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i])*(r_g_vector[j] - state_old -> entropy_densities[j*NO_OF_SCALARS_H + i]/state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]);
			beta[j] = state_old -> temperature_gas[j*NO_OF_SCALARS_H + i]/(c_g_v*state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]);
			// partial derivatives of s = s(rho, stilde)
			gamma[j] = -state_old -> entropy_densities[j*NO_OF_SCALARS_H + i]/pow(state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i], 2);
			delta[j] = 1/state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			// dividing by the volume:
			alpha[j] = alpha[j]/grid -> volume[i + j*NO_OF_SCALARS_H];
			beta[j] = beta[j]/grid -> volume[i + j*NO_OF_SCALARS_H];
			gamma[j] = gamma[j]/grid -> volume[i + j*NO_OF_SCALARS_H];
			delta[j] = delta[j]/grid -> volume[i + j*NO_OF_SCALARS_H];
		}
		
		// determining the upper and lower weights as well as the interface values
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
            lower_index = i + (j + 1)*NO_OF_SCALARS_H;
            upper_index = i + j*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            a[j] = upper_volume/total_volume;
            b[j] = lower_volume/total_volume;
			temp_old_interface_values[j] = a[j]*state_old -> temperature_gas[upper_index]
			+ b[j]*state_old -> temperature_gas[lower_index];
			spec_entropy_interface_new[j]
			= a[j]*spec_entropy_new[j] + b[j]*spec_entropy_new[j + 1];
			density_interface_old[j]
			= a[j]*state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			+ b[j]*state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (j + 1)*NO_OF_SCALARS_H + i];
			density_interface_explicit[j]
			= a[j]*(state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			+ delta_t*state_tendency -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i])
			+ b[j]*(state_old -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (j + 1)*NO_OF_SCALARS_H + i]
			+ delta_t*state_tendency -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (j + 1)*NO_OF_SCALARS_H + i]);
		}
		
		// filling up the coefficient vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{			
			d_vector[j] = -c_g_p*alpha[j] - c_g_p*(beta[j] + beta[j + 1])*spec_entropy_interface_new[j] - c_g_p*alpha[j + 1]
			- (grid -> z_scalar[i + j*NO_OF_SCALARS_H] - grid -> z_scalar[i + (j + 1)*NO_OF_SCALARS_H])/(impl_pgrad_weight*pow(delta_t, 2))*2/grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
			// right hand side
			r_vector[j] = -(grid -> z_scalar[i + j*NO_OF_SCALARS_H] - grid -> z_scalar[i + (j + 1)*NO_OF_SCALARS_H])/(impl_pgrad_weight*pow(delta_t, 2))
			*(-density_interface_explicit[j]*state_old -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER])/(density_interface_old[j])
			- (state_old -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] + delta_t*state_tendency -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER])
			*(grid -> z_scalar[i + j*NO_OF_SCALARS_H] - grid -> z_scalar[i + (j + 1)*NO_OF_SCALARS_H])/(impl_pgrad_weight*pow(delta_t, 2))
			+ c_g_p/delta_t*(diagnostics -> temperature_gas_explicit[i + j*NO_OF_SCALARS_H] - diagnostics -> temperature_gas_explicit[i + (j + 1)*NO_OF_SCALARS_H])
			- temp_old_interface_values[j]/delta_t*(spec_entropy_explicit[j] - spec_entropy_explicit[j + 1]);
		}
		for (j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			c_vector[j] = c_g_p*alpha[j + 1] + c_g_p*beta[j + 1]*spec_entropy_interface_new[j];
			e_vector[j] = c_g_p*alpha[j] + c_g_p*beta[j]*spec_entropy_interface_new[j];
		}
		
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_LAYERS - 1);
		
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
			solution_vector[j] = solution_vector[j]/(1 + delta_t*damping_coeff);
		}
		
		/*
		Writing the result into the new state.
		--------------------------------------
		*/
		// mass density
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				= density_explicit[j] + delta_t*(solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				= density_explicit[j] + delta_t*(-solution_vector[j - 1])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else
			{
				state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				= density_explicit[j] + delta_t*(-solution_vector[j - 1] + solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
		}
		// entropy density
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				state_new -> entropy_densities[j*NO_OF_SCALARS_H + i]
				= entropy_density_explicit[j] + delta_t*(entropy_density_explicit[j]*solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				state_new -> entropy_densities[j*NO_OF_SCALARS_H + i]
				= entropy_density_explicit[j] + delta_t*(-entropy_density_explicit[j - 1]*solution_vector[j - 1])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
			else
			{
				state_new -> entropy_densities[j*NO_OF_SCALARS_H + i]
				= entropy_density_explicit[j] + delta_t*(-entropy_density_explicit[j - 1]*solution_vector[j - 1] + entropy_density_explicit[j]*solution_vector[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
			}
		}
		// vertical velocity
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			density_interface_new
			= a[j]*state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
			+ b[j]*state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (j + 1)*NO_OF_SCALARS_H + i];
			state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]
			= (2*solution_vector[j]/grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - density_interface_new*state_old -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i])
			/density_interface_old[j];
		}
		// temperature
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> temperature_gas[j*NO_OF_SCALARS_H + i]
			= diagnostics -> temperature_gas_explicit[j*NO_OF_SCALARS_H + i]
			+ grid -> volume[i + j*NO_OF_SCALARS_H]*alpha[j]*(state_new -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] - density_explicit[j])
			+ grid -> volume[i + j*NO_OF_SCALARS_H]*beta[j]*(state_new -> entropy_densities[j*NO_OF_SCALARS_H + i] - entropy_density_explicit[j]);
		}
	}
	return 0;
}

int three_band_solver_gen_densitites(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, Grid *grid)
{
	// Vertical advection of generalized densities (of tracers) with 3-band matrices.
	// mass densities, entropy densities, density x temperatures
	int no_of_relevant_constituents, constituent_index_offset;
	double impl_weight, expl_weight;
	impl_weight = 0.5;
	expl_weight = 1 - impl_weight;
	for (int quantity_id = 0; quantity_id < 2; ++quantity_id)
	{
		no_of_relevant_constituents = 0;
		constituent_index_offset = 0;
		// mass densities
		if (quantity_id == 0)
		{
			// all constituents have a mass density
			no_of_relevant_constituents = NO_OF_CONSTITUENTS - 1;
			constituent_index_offset = 0;
		}
		// density x temperature fields
		if (quantity_id == 1)
		{
			constituent_index_offset = 0;
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
		for (int k = constituent_index_offset; k < constituent_index_offset + no_of_relevant_constituents; ++k)
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
					double a[NO_OF_LAYERS - 1];
					double b[NO_OF_LAYERS - 1];
					double solution_vector[NO_OF_LAYERS];
					double density_gas_value, density_gas_value_old;
					double density_old_at_interface, area, upper_volume, lower_volume, total_volume;
					int j, lower_index, upper_index;
					
					// diagnozing the vertical fluxes
					for (j = 0; j < NO_OF_LAYERS - 1; ++j)
					{
						vertical_flux_vector_impl[j] = state_old -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
						vertical_flux_vector_rhs[j] = state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
						// preparing the vertical interpolation
						lower_index = i + (j + 1)*NO_OF_SCALARS_H;
						upper_index = i + j*NO_OF_SCALARS_H;
						upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
						lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
						total_volume = upper_volume + lower_volume;
						a[j] = upper_volume/total_volume;
						b[j] = lower_volume/total_volume;
						// For condensed constituents, a sink velocity must be added.
						if (k < NO_OF_CONDENSED_CONSTITUENTS)
						{
							// determining the density of the gas at the interface
							density_gas_value = a[j]*density_gas(state_new, j*NO_OF_SCALARS_H + i) + b[j]*density_gas(state_new, (j + 1)*NO_OF_SCALARS_H + i);
							density_gas_value_old = a[j]*density_gas(state_old, j*NO_OF_SCALARS_H + i) + b[j]*density_gas(state_old, (j + 1)*NO_OF_SCALARS_H + i);
							vertical_flux_vector_impl[j] -= ret_sink_velocity(0, 0.001, density_gas_value_old);
							vertical_flux_vector_rhs[j] -= ret_sink_velocity(0, 0.001, density_gas_value);
						}
						// multiplying the vertical velocity by the area
						area = grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
						vertical_flux_vector_impl[j] = area*vertical_flux_vector_impl[j];
						vertical_flux_vector_rhs[j] = area*vertical_flux_vector_rhs[j];
						// old density at the interface
						density_old_at_interface
						= a[j]*state_old -> mass_densities[k*NO_OF_SCALARS + upper_index]
						+ b[j]*state_old -> mass_densities[k*NO_OF_SCALARS + lower_index];
						vertical_flux_vector_rhs[j] = density_old_at_interface*vertical_flux_vector_rhs[j];
					}
					
					/*
					Now we proceed to solve the vertical tridiagonal problems.
					*/
					// filling up the original vectors
					for (j = 0; j < NO_OF_LAYERS - 1; ++j)
					{
						c_vector[j] = impl_weight*a[j]*delta_t/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
						e_vector[j] = -impl_weight*b[j]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
					}
					for (j = 0; j < NO_OF_LAYERS; ++j)
					{
						if (j == 0)
						{
							d_vector[j] = 1
							- impl_weight*a[j]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[0];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							d_vector[j] = 1
							+ impl_weight*b[j - 1]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j - 1];
						}
						else
						{
							d_vector[j] = 1
							+ impl_weight*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]
							*(b[j - 1]*vertical_flux_vector_impl[j - 1] - a[j]*vertical_flux_vector_impl[j]);
						}
						// the explicit component
						// mass densities
						if (quantity_id == 0)
						{
							r_vector[j] =
							state_old -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
							+ delta_t*state_tendency -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
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
					thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_LAYERS);
					for (j = 0; j < NO_OF_LAYERS; ++j)
					{
						// limiter: none of the densities may be negative
						if (solution_vector[j] < 0)
						{
							solution_vector[j] = 0;
						}
						if (quantity_id == 0)
						{
							state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
						}
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
	// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	double c_prime_vector[solution_length - 1];
	double d_prime_vector[solution_length];
	if (d_vector[0] != 0)
	{
		c_prime_vector[0] = e_vector[0]/d_vector[0];
	}
	else
	{
		c_prime_vector[0] = 0;
	}
	for (int j = 1; j < solution_length - 1; ++j)
	{
		if (d_vector[j] - c_prime_vector[j - 1]*c_vector[j - 1] != 0)
		{
			c_prime_vector[j] = e_vector[j]/(d_vector[j] - c_prime_vector[j - 1]*c_vector[j - 1]);
		}
		else
		{
			c_prime_vector[j] = 0;
		}
	}
	if (d_vector[0] != 0)
	{
		d_prime_vector[0] = r_vector[0]/d_vector[0];
	}
	else
	{
		d_prime_vector[0] = 0;
	}
	for (int j = 1; j < solution_length; ++j)
	{
		if (d_vector[j] - c_prime_vector[j - 1]*c_vector[j - 1] != 0)
		{
			d_prime_vector[j] = (r_vector[j] - d_prime_vector[j - 1]*c_vector[j - 1])/(d_vector[j] - c_prime_vector[j - 1]*c_vector[j - 1]);
		}
		else
		{
			d_prime_vector[j] = 0;
		}
	}
	solution_vector[solution_length - 1] = d_prime_vector[solution_length - 1];
	for (int j = solution_length - 2; j >= 0; --j)
	{
		solution_vector[j] = d_prime_vector[j] - c_prime_vector[j]*solution_vector[j + 1];
	}
	return 0;
}



