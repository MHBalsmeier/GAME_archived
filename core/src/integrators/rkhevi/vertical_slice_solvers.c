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

int three_band_solver_ver_sound_waves(State *state_old, State *state_tendency, State *state_new, Diagnostics *diagnostics, Config_info *config_info, double delta_t, Grid *grid)
{
	double delta_z, upper_volume, lower_volume, total_volume, damping_coeff, damping_coeff_max, damping_start_height, z_above_damping, damping_start_height_over_toa;
	// This is for Klemp (2008).
	get_damping_layer_properties(&damping_start_height_over_toa, &damping_coeff_max);
	damping_start_height = damping_start_height_over_toa*grid -> z_vector[0];
	int upper_index, lower_index, j;
	double impl_pgrad_weight = get_impl_thermo_weight();
	#pragma omp parallel for private(upper_index, lower_index, j, delta_z, upper_volume, lower_volume, total_volume, damping_coeff)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[2*NO_OF_LAYERS - 2];
		double b_vector[2*NO_OF_LAYERS - 1];
		double c_vector[2*NO_OF_LAYERS - 2];
		double d_vector[2*NO_OF_LAYERS - 1];
		double solution_vector[2*NO_OF_LAYERS - 1];
		double upper_weights_vector[NO_OF_LAYERS - 1];
		double lower_weights_vector[NO_OF_LAYERS - 1];
		double temp_interface_values[NO_OF_LAYERS - 1];
		double c_g_v_vector[NO_OF_LAYERS];
		double c_g_p_vector[NO_OF_LAYERS];
		double r_g_vector[NO_OF_LAYERS];
		double c_g_p_interface_values[NO_OF_LAYERS - 1];
		// determining the properties of the gas phase in the grid boxes
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			c_g_v_vector[j] = spec_heat_cap_diagnostics_v(state_new, i + j*NO_OF_SCALARS_H, config_info);
			c_g_p_vector[j] = spec_heat_cap_diagnostics_p(state_new, i + j*NO_OF_SCALARS_H, config_info);
			r_g_vector[j] = gas_constant_diagnostics(state_new, i + j*NO_OF_SCALARS_H, config_info);
		}
		// determining the upper and lower weights as well as the c_g_v interface values
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
            lower_index = i + (j + 1)*NO_OF_SCALARS_H;
            upper_index = i + j*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            upper_weights_vector[j] = upper_volume/total_volume;
            lower_weights_vector[j] = lower_volume/total_volume;
			temp_interface_values[j] = upper_weights_vector[j]*state_new -> temperature_gas[upper_index]
			+ lower_weights_vector[j]*state_new -> temperature_gas[lower_index];
			c_g_p_interface_values[j] = upper_weights_vector[j]*c_g_p_vector[j] +
			lower_weights_vector[j]*c_g_p_vector[j + 1];
		}
		// filling up the original vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{			
			// determining the elements of a_vector
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			a_vector[2*j] = delta_t*impl_pgrad_weight*c_g_p_interface_values[j]/delta_z;
			a_vector[2*j + 1] = delta_t*((r_g_vector[j]/c_g_v_vector[j] - 1)*state_new -> temperature_gas[i + (j + 1)*NO_OF_SCALARS_H] + temp_interface_values[j])*
			grid -> area[i + j*NO_OF_VECTORS_PER_LAYER]/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H];
			
			// determining the elements of c_vector
			c_vector[2*j] = -delta_t*((r_g_vector[j]/c_g_v_vector[j] - 1)*state_new -> temperature_gas[i + j*NO_OF_SCALARS_H] + temp_interface_values[j])*
			grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]/grid -> volume[i + j*NO_OF_SCALARS_H];
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			c_vector[2*j + 1] = -delta_t*impl_pgrad_weight*c_g_p_interface_values[j]/delta_z;
		}
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			b_vector[2*j] = 1;
			// Klemp (2008) upper boundary layer
			z_above_damping = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - damping_start_height;
			if (z_above_damping < 0)
			{
				damping_coeff = 0;
			}
			else
			{
				damping_coeff = damping_coeff_max*pow(sin(0.5*M_PI*z_above_damping/(grid -> z_vector[0] - damping_start_height)), 2);
			}
			b_vector[2*j + 1] = 1 + delta_t*damping_coeff;
			d_vector[2*j] = diagnostics -> temperature_gas_explicit[j*NO_OF_SCALARS_H + i];
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			d_vector[2*j + 1] =
			state_old -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]
			// explicit tendency
			+ delta_t*state_tendency -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
		}
		b_vector[2*NO_OF_LAYERS - 2] =  1;
		d_vector[2*NO_OF_LAYERS - 2] = diagnostics -> temperature_gas_explicit[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, solution_vector, 2*NO_OF_LAYERS - 1);
		// writing the result into the new state
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			state_new -> temperature_gas[j*NO_OF_SCALARS_H + i] = solution_vector[2*j];
			state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[2*j + 1];
		}
		state_new -> temperature_gas[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i] = solution_vector[2*(NO_OF_LAYERS - 1)];
	}
	return 0;
}

int three_band_solver_gen_densitites(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, Grid *grid)
{
	// Vertical constituent advection with 3-band matrices.
	// procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	// mass densities, entropy densities, density x temperatures
	int no_of_relevant_constituents, constituent_index_offset;
	for (int quantity_id = 0; quantity_id < 3; ++quantity_id)
	{
		no_of_relevant_constituents = 0;
		constituent_index_offset = 0;
		// mass densities
		if (quantity_id == 0)
		{
			// all constituents have a mass density
			no_of_relevant_constituents = NO_OF_CONSTITUENTS;
			constituent_index_offset = 0;
		}
		// entropy densities
		if (quantity_id == 1)
		{
			// in this case, all gaseous constituents have an entropy density
			if (config_info -> assume_lte == 0)
			{
				constituent_index_offset = NO_OF_CONDENSED_CONSTITUENTS;
				no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
			}
			// in this case, only one constituent has an entropy density
			if (config_info -> assume_lte == 1)
			{
				constituent_index_offset = NO_OF_CONDENSED_CONSTITUENTS;
				no_of_relevant_constituents = 1;
			}
		}
		// density x temperature fields
		if (quantity_id == 2)
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
			// loop over all columns
			#pragma omp parallel for
			for (int i = 0; i < NO_OF_SCALARS_H; ++i)
			{
				// for meanings of these vectors look into the definition of the function thomas_algorithm
				double a_vector[NO_OF_LAYERS - 1];
				double b_vector[NO_OF_LAYERS];
				double c_vector[NO_OF_LAYERS - 1];
				double d_vector[NO_OF_LAYERS];
				double vertical_flux_vector_impl[NO_OF_LAYERS - 1];
				double vertical_flux_vector_rhs[NO_OF_LAYERS - 1];
				double upper_weights_vector[NO_OF_LAYERS - 1];
				double lower_weights_vector[NO_OF_LAYERS - 1];
				double solution_vector[NO_OF_LAYERS];
				double density_gas_value, density_gas_value_old;
				double density_new_at_interface, density_old_at_interface, area, upper_volume, lower_volume, total_volume;
				int j, lower_index, upper_index;
				// diagnozing the vertical flux
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
					upper_weights_vector[j] = upper_volume/total_volume;
					lower_weights_vector[j] = lower_volume/total_volume;
					// For condensed constituents, a sink velocity must be added.
					if (k < NO_OF_CONDENSED_CONSTITUENTS)
					{
						// determining the density of the gas at the interface
						density_gas_value = upper_weights_vector[j]*density_gas(state_new, j*NO_OF_SCALARS_H + i) + lower_weights_vector[j]*density_gas(state_new, (j + 1)*NO_OF_SCALARS_H + i);
						density_gas_value_old = upper_weights_vector[j]*density_gas(state_old, j*NO_OF_SCALARS_H + i) + lower_weights_vector[j]*density_gas(state_old, (j + 1)*NO_OF_SCALARS_H + i);
						vertical_flux_vector_impl[j] -= ret_sink_velocity(0, 0.001, density_gas_value_old);
						vertical_flux_vector_rhs[j] -= ret_sink_velocity(0, 0.001, density_gas_value);
					}
					// multiplying the vertical velocity by the area
					area = grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
					vertical_flux_vector_impl[j] = area*vertical_flux_vector_impl[j];
					vertical_flux_vector_rhs[j] = area*vertical_flux_vector_rhs[j];
					// old density at the interface
					density_old_at_interface = upper_weights_vector[j]*state_old -> mass_densities[k*NO_OF_SCALARS + upper_index]
					+ lower_weights_vector[j]*state_old -> mass_densities[k*NO_OF_SCALARS + lower_index];
					// for entropy densities, the vertical_flux_density is the mass flux density
					if (quantity_id == 1)
					{
						// new density at the interface
						density_new_at_interface =
						upper_weights_vector[j]*state_new -> mass_densities[k*NO_OF_SCALARS + upper_index]
						+ lower_weights_vector[j]*state_new -> mass_densities[k*NO_OF_SCALARS + lower_index];
						vertical_flux_vector_impl[j] = density_old_at_interface*vertical_flux_vector_impl[j];
						vertical_flux_vector_rhs[j] =
						// the old specific entropy at the interface
						(upper_weights_vector[j]*state_old -> entropy_densities[(k - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + upper_index]/state_old -> mass_densities[k*NO_OF_SCALARS + upper_index]
						+ lower_weights_vector[j]*state_old -> entropy_densities[(k - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + lower_index]/state_old -> mass_densities[k*NO_OF_SCALARS + lower_index])
						 // the mass flux density used in ther vertical mass flux divergence
						*density_new_at_interface*vertical_flux_vector_rhs[j];
					}
					else
					{
						vertical_flux_vector_rhs[j] = density_old_at_interface*vertical_flux_vector_rhs[j];
					}
				}
				// filling up the original vectors
				for (j = 0; j < NO_OF_LAYERS - 1; ++j)
				{
					// entropy densities
					if (quantity_id == 1)
					{
						a_vector[j] = 0.5*upper_weights_vector[j]*delta_t/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
						c_vector[j] = -0.5*lower_weights_vector[j]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
					}
					else
					{
						a_vector[j] = 0.5*upper_weights_vector[j]*delta_t/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
						c_vector[j] = -0.5*lower_weights_vector[j]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
					}
				}
				for (j = 0; j < NO_OF_LAYERS; ++j)
				{
					// entropy densities
					if (quantity_id == 1)
					{
						if (j == 0)
						{
							b_vector[j] = state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] - 0.5*upper_weights_vector[j]*delta_t
							/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[0];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							b_vector[j] = state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + 0.5*lower_weights_vector[j - 1]*delta_t
							/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j - 1];
						}
						else
						{
							b_vector[j] = state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + 0.5*delta_t
							/grid -> volume[i + j*NO_OF_SCALARS_H]*(lower_weights_vector[j - 1]*vertical_flux_vector_impl[j - 1] - upper_weights_vector[j]*vertical_flux_vector_impl[j]);
						}
					}
					else
					{
						if (j == 0)
						{
							b_vector[j] = 1 - 0.5*upper_weights_vector[j]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[0];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							b_vector[j] = 1 + 0.5*lower_weights_vector[j - 1]*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j - 1];
						}
						else
						{
							b_vector[j] = 1 + 0.5*delta_t/grid -> volume[i + j*NO_OF_SCALARS_H]
							*(lower_weights_vector[j - 1]*vertical_flux_vector_impl[j - 1] - upper_weights_vector[j]*vertical_flux_vector_impl[j]);
						}
					}
					// the explicit component
					// mass densities
					if (quantity_id == 0)
					{
						d_vector[j] =
						state_old -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
						+ delta_t*state_tendency -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
					}
					// entropy densities
					if (quantity_id == 1)
					{
						d_vector[j] =
						state_old -> entropy_densities[(k - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
						+ delta_t*state_tendency -> entropy_densities[(k - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
					}
					// density x temperatues
					if (quantity_id == 2)
					{
						d_vector[j] =
						state_old -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
						+ delta_t*state_tendency -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
					}
					// adding the explicit part of the vertical flux divergence
					if (j == 0)
					{
						d_vector[j] += 0.5*delta_t*vertical_flux_vector_rhs[j]/grid -> volume[j*NO_OF_SCALARS_H + i];
					}
					else if (j == NO_OF_LAYERS - 1)
					{
						d_vector[j] += -0.5*delta_t*vertical_flux_vector_rhs[j - 1]/grid -> volume[j*NO_OF_SCALARS_H + i];
					}
					else
					{
						d_vector[j] += 0.5*delta_t*(-vertical_flux_vector_rhs[j - 1] + vertical_flux_vector_rhs[j])/grid -> volume[j*NO_OF_SCALARS_H + i];
					}
				}
				thomas_algorithm(a_vector, b_vector, c_vector, d_vector, solution_vector, NO_OF_LAYERS);
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
					 	// here, the solution vector actually contains the specific entropy, that's why it needs to be multiplied by the mass density
						state_new -> entropy_densities[(k - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] =
						state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]*solution_vector[j];
					}
					if (quantity_id == 2)
					{
						state_new -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
					}
				}
			}
		}
	}
	return 0;
}

int thomas_algorithm(double a_vector[], double b_vector[], double c_vector[], double d_vector[], double solution_vector[], int solution_length)
{
	// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	double c_prime_vector[solution_length - 1];
	double d_prime_vector[solution_length];
	if (b_vector[0] != 0)
	{
		c_prime_vector[0] = c_vector[0]/b_vector[0];
	}
	else
	{
		c_prime_vector[0] = 0;
	}
	for (int j = 1; j < solution_length - 1; ++j)
	{
		if (b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1] != 0)
		{
			c_prime_vector[j] = c_vector[j]/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1]);
		}
		else
		{
			c_prime_vector[j] = 0;
		}
	}
	if (b_vector[0] != 0)
	{
		d_prime_vector[0] = d_vector[0]/b_vector[0];
	}
	else
	{
		d_prime_vector[0] = 0;
	}
	for (int j = 1; j < solution_length; ++j)
	{
		if (b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1] != 0)
		{
			d_prime_vector[j] = (d_vector[j] - d_prime_vector[j - 1]*a_vector[j - 1])/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1]);
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



