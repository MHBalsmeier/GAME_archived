/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
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

int thomas_algorithm(double [], double [], double [], double [], double [], double [], double [], int);
int lu_5band_solver(double [], double [], double [], double [], double [], double [], double [], int);
int sign(double);

int three_band_solver_ver_hor_vel_adv(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Semi-implicit vertical advection of vertical momentum (Crank-Nicolson with variable implicit weight).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double delta_z;
	int i, j;
	double impl_u_vadv_weight = get_impl_u_vadv_weight();
	double expl_u_vadv_weight = 1 - impl_u_vadv_weight;
	#pragma omp parallel for private(delta_z, j)
	for (i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[NO_OF_LAYERS - 1];
		double b_vector[NO_OF_LAYERS];
		double c_vector[NO_OF_LAYERS - 1];
		double d_vector[NO_OF_LAYERS];
		double c_prime_vector[NO_OF_LAYERS - 1];
		double d_prime_vector[NO_OF_LAYERS];
		double vertical_velocity[NO_OF_LAYERS];
		double solution_vector[NO_OF_LAYERS];
		// diagnozing the vertical velocity
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			recov_hor_ver_pri(state_new -> velocity_gas, j, i, &vertical_velocity[j], grid);
		}
		// filling up b vector and d vector
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 + 0.5*delta_t*vertical_velocity[j]/(2*delta_z);
				// right hand side
				// old value
				d_vector[j] = state_old -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// the old time step part of Crank-Nicolson
				- expl_u_vadv_weight*delta_t*vertical_velocity[j]
				*(state_new -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - state_new -> velocity_gas[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i])
				/delta_z;
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				b_vector[j] = 1;
				delta_z = 2*(grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]);
				// right hand side
				// old value
				d_vector[j] = state_old -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// the old time step part of Crank-Nicolson
				- expl_u_vadv_weight*delta_t*vertical_velocity[j]*state_new -> velocity_gas[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i]/delta_z;
			}
			else
			{
				b_vector[j] = 1;
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
				// right hand side
				// old value
				d_vector[j] = state_old -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// the old time step part of Crank-Nicolson
				- expl_u_vadv_weight*delta_t*vertical_velocity[j]
				*(state_new -> velocity_gas[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - state_new -> velocity_gas[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i])/delta_z;
			}
		}
		// filling up a vector and c vector
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			if (j == NO_OF_LAYERS - 2)
			{
				delta_z = 2*(grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i]);
			}
			else
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			}
			a_vector[j] = impl_u_vadv_weight*vertical_velocity[j + 1]*delta_t/(2*delta_z);
			if (j == 0)
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			}
			else
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			}
			c_vector[j] = -impl_u_vadv_weight*vertical_velocity[j]*delta_t/(2*delta_z);
		}
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		// writing the solution into the new state
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[j];
		}
	}
	return 0;
}

int three_band_solver_ver_sound_waves(State *state_new, Diagnostics *diagnostics, Interpolation_info *interpolation_info, double delta_t, Grid *grid)
{
	double delta_z, upper_volume, lower_volume, total_volume, damping_coeff, damping_coeff_max, damping_start_height, z_above_damping, damping_start_height_over_toa;
	// This is for Klemp (2008).
	get_damping_layer_properties(&damping_start_height_over_toa, &damping_coeff_max);
	damping_start_height = damping_start_height_over_toa*grid -> z_vector[0];
	int upper_index, lower_index, j, gaseous_constituent_id;
	double impl_pgrad_weight = 1 - get_expl_pgrad_weight();
	double impl_vert_compression_weight = get_impl_vert_compression_weight();
	#pragma omp parallel for private(upper_index, lower_index, delta_z, upper_volume, lower_volume, total_volume, damping_coeff, z_above_damping, j, gaseous_constituent_id)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function lu_5band_solver
		double a_vector[2*NO_OF_LAYERS - 2];
		double b_vector[2*NO_OF_LAYERS - 1];
		double c_vector[2*NO_OF_LAYERS - 2];
		double d_vector[2*NO_OF_LAYERS - 1];
		double solution_vector[2*NO_OF_LAYERS - 1];
		double upper_weights_vector[NO_OF_LAYERS - 1];
		double lower_weights_vector[NO_OF_LAYERS - 1];
		double vertical_entropy_vector[NO_OF_LAYERS];
		double vertical_entropy_gradient_component[NO_OF_LAYERS - 1];
		double vertical_entropy_gradient[NO_OF_GASEOUS_CONSTITUENTS*(NO_OF_LAYERS - 1)];
		double temp_interface_values[NO_OF_LAYERS - 1];
		double c_g_v_vector[NO_OF_LAYERS];
		double c_g_p_vector[NO_OF_LAYERS];
		double r_g_vector[NO_OF_LAYERS];
		double c_g_p_interface_values[NO_OF_LAYERS - 1];
		double density_gas_vector[NO_OF_LAYERS];
		double vertical_velocity_divergence[NO_OF_LAYERS];
		double l_vector[2*NO_OF_LAYERS - 3];
		double u_vector[2*NO_OF_LAYERS - 3];
		// determining the properties of the gas phase in the grid boxes and the vertical velocity divergence
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			c_g_v_vector[j] = spec_heat_cap_diagnostics_v(state_new, i + j*NO_OF_SCALARS_H);
			c_g_p_vector[j] = spec_heat_cap_diagnostics_p(state_new, i + j*NO_OF_SCALARS_H);
			r_g_vector[j] = gas_constant_diagnostics(state_new, i + j*NO_OF_SCALARS_H);	
			density_gas_vector[j] = density_gas(state_new, i + j*NO_OF_SCALARS_H);
			if (j == 0)
			{
				vertical_velocity_divergence[j] =
				-grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]*interpolation_info -> velocity_gas_prior_rk[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]
				/grid -> volume[i + j*NO_OF_SCALARS_H];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				vertical_velocity_divergence[j] = 
				grid -> area[i + j*NO_OF_VECTORS_PER_LAYER]*interpolation_info -> velocity_gas_prior_rk[i + j*NO_OF_VECTORS_PER_LAYER]
				/grid -> volume[i + j*NO_OF_SCALARS_H];
			}
			else
			{
				vertical_velocity_divergence[j] = (
				grid -> area[i + j*NO_OF_VECTORS_PER_LAYER]*interpolation_info -> velocity_gas_prior_rk[i + j*NO_OF_VECTORS_PER_LAYER]
				- grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]*interpolation_info -> velocity_gas_prior_rk[i + (j + 1)*NO_OF_VECTORS_PER_LAYER])
				/grid -> volume[i + j*NO_OF_SCALARS_H];
			}
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
			// initializing vertical_entropy_gradient with zeros
			for (gaseous_constituent_id = 0; gaseous_constituent_id < NO_OF_GASEOUS_CONSTITUENTS; ++gaseous_constituent_id)
			{
				vertical_entropy_gradient[gaseous_constituent_id*(NO_OF_LAYERS - 1) + j] = 0;
			}
		}
		// determining the specific entropy gradient
		// loop over all gaseous constituents
		for (gaseous_constituent_id = 0; gaseous_constituent_id < NO_OF_GASEOUS_CONSTITUENTS; ++gaseous_constituent_id)
		{
			// determining the specific entropy of the constituent at hand
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{	
				if (state_new -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + gaseous_constituent_id)*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H] != 0)
				{
					vertical_entropy_vector[j] = state_new -> entropy_densities[(NO_OF_CONDENSED_CONSTITUENTS + gaseous_constituent_id)*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H]/
					state_new -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + gaseous_constituent_id)*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H];
				}
				else
				{
					vertical_entropy_vector[j] = 0;
				}
			}
			// taking the gradient
			grad_v_scalar_column(vertical_entropy_vector, vertical_entropy_gradient_component, i, grid);
			// adding to the density weighted sum of all specific entropy gradients
			for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				vertical_entropy_gradient[gaseous_constituent_id*(NO_OF_LAYERS - 1) + j] = vertical_entropy_gradient_component[j];
			}
		}
		// filling up the original vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{			
			// determining the elements of a_vector
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			a_vector[2*j] = delta_t*impl_pgrad_weight*c_g_p_interface_values[j]/delta_z;
			for (gaseous_constituent_id = 0; gaseous_constituent_id < NO_OF_GASEOUS_CONSTITUENTS; ++gaseous_constituent_id)
			{
				a_vector[2*j] +=
				- delta_t*upper_weights_vector[j]*impl_pgrad_weight
				*state_new -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + gaseous_constituent_id)*NO_OF_SCALARS + i + j*NO_OF_SCALARS_H]/density_gas_vector[j]
				*vertical_entropy_gradient[gaseous_constituent_id*(NO_OF_LAYERS - 1) + j];
			}
			a_vector[2*j + 1] = delta_t*((impl_vert_compression_weight*r_g_vector[j]/c_g_v_vector[j] - 1)*state_new -> temperature_gas[i + (j + 1)*NO_OF_SCALARS_H] + temp_interface_values[j])*
			grid -> area[i + j*NO_OF_VECTORS_PER_LAYER]/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H];
			
			// determining the elements of c_vector
			c_vector[2*j] = -delta_t*((impl_vert_compression_weight*r_g_vector[j]/c_g_v_vector[j] - 1)*state_new -> temperature_gas[i + j*NO_OF_SCALARS_H] + temp_interface_values[j])*
			grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]/grid -> volume[i + j*NO_OF_SCALARS_H];
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			c_vector[2*j + 1] = -delta_t*impl_pgrad_weight*c_g_p_interface_values[j]/delta_z;
			for (gaseous_constituent_id = 0; gaseous_constituent_id < NO_OF_GASEOUS_CONSTITUENTS; ++gaseous_constituent_id)
			{
				c_vector[2*j + 1] +=
				- delta_t*lower_weights_vector[j]*impl_pgrad_weight
				*state_new -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + gaseous_constituent_id)*NO_OF_SCALARS + i + (j + 1)*NO_OF_SCALARS_H]/density_gas_vector[j + 1]
				*vertical_entropy_gradient[gaseous_constituent_id*(NO_OF_LAYERS - 1) + j];
			}
		}
		for (int j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			// determining the elements of l_vector
			l_vector[2*j] = 0;
			l_vector[2*j + 1] = 0;
			// determining the elements of u_vector
			u_vector[2*j] = 0;
			u_vector[2*j + 1] = 0;
		}
		l_vector[2*NO_OF_LAYERS - 4] = 0;
		u_vector[2*NO_OF_LAYERS - 4] = 0;
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			b_vector[2*j] = 1 + delta_t*(1 - impl_vert_compression_weight)*r_g_vector[j]/c_g_v_vector[j]*vertical_velocity_divergence[j];
			b_vector[2*j + 1] = 1;
			d_vector[2*j] = diagnostics -> temperature_gas_explicit[j*NO_OF_SCALARS_H + i];
			d_vector[2*j + 1] = state_new -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
		}
		b_vector[2*NO_OF_LAYERS - 2] =  1 + delta_t*(1 - impl_vert_compression_weight)*r_g_vector[NO_OF_LAYERS - 1]/c_g_v_vector[NO_OF_LAYERS - 1]*vertical_velocity_divergence[NO_OF_LAYERS - 1];
		d_vector[2*NO_OF_LAYERS - 2] = diagnostics -> temperature_gas_explicit[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
		// calling the algorithm to solve the system of linear equations
		lu_5band_solver(a_vector, b_vector, c_vector, d_vector, l_vector, u_vector, solution_vector, 2*NO_OF_LAYERS - 1);
		// writing the result into the new state
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
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
			state_new -> temperature_gas[j*NO_OF_SCALARS_H + i] = solution_vector[2*j];
			state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[2*j + 1]/(1 + delta_t*damping_coeff);
			interpolation_info -> velocity_gas_prior_rk[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
		}
		state_new -> temperature_gas[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i] = solution_vector[2*(NO_OF_LAYERS - 1)];
	}
	return 0;
}

int three_band_solver_gen_densitites(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	// Vertical constituent advection with 3-band matrices (Euler implicit).
	// procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	// for meanings of these vectors look into the definition of the function thomas_algorithm
	double a_vector[NO_OF_LAYERS - 1];
	double b_vector[NO_OF_LAYERS];
	double c_vector[NO_OF_LAYERS - 1];
	double d_vector[NO_OF_LAYERS];
	double c_prime_vector[NO_OF_LAYERS - 1];
	double d_prime_vector[NO_OF_LAYERS];
	double vertical_flux_vector[NO_OF_LAYERS - 1];
	double solution_vector[NO_OF_LAYERS];
	double density_gas_value;
	int no_of_relevant_constituents = 0;
	double new_density_at_interface, area, upper_volume, lower_volume, total_volume, upper_weight, lower_weight;
	int j, lower_index, upper_index;
	// mass densities, entropy densities, density x temperatures
	for (int quantity_id = 0; quantity_id < 3; ++quantity_id)
	{
		// all constituents have a mass density
		if (quantity_id == 0)
		{
			no_of_relevant_constituents = NO_OF_CONSTITUENTS;
		}
		// all constituents have an entropy density
		if (quantity_id == 1)
		{
			no_of_relevant_constituents = NO_OF_CONSTITUENTS;
		}
		// only the condensed constituents have a density x temperature field
		if (quantity_id == 2)
		{
			no_of_relevant_constituents = NO_OF_CONDENSED_CONSTITUENTS;
		}
		// loop over all relevant constituents
		for (int k = 0; k < no_of_relevant_constituents; ++k)
		{
			#pragma omp parallel for private(area, j, new_density_at_interface, density_gas_value, lower_index, upper_index, upper_volume, lower_volume, total_volume, upper_weight, lower_weight)
			for (int i = 0; i < NO_OF_SCALARS_H; ++i)
			{
				// diagnozing the vertical flux
				for (j = 0; j < NO_OF_LAYERS - 1; ++j)
				{
					vertical_flux_vector[j] = state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
					// For condensed contituents, a sink velocity must be added.
					if (k < NO_OF_CONDENSED_CONSTITUENTS)
					{
						// determining the density of the gas
						density_gas_value = 0.5*(density_gas(state_new, j*NO_OF_SCALARS_H + i) + density_gas(state_new, (j + 1)*NO_OF_SCALARS_H + i));
						
						// The solid case.
						if (k < NO_OF_SOLID_CONSTITUENTS)
						{
							vertical_flux_vector[j] -= ret_sink_velocity(0, 0.001, density_gas_value);
						}
						// The liquid case.
						else
						{
							vertical_flux_vector[j] -= ret_sink_velocity(1, 0.001, density_gas_value);
						}
					}
					area = grid -> area[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
					vertical_flux_vector[j] = area*vertical_flux_vector[j];
					// for entropy densities, the vertical flux vector is the mass flux density
					if (quantity_id == 1)
					{
						lower_index = i + (j + 1)*NO_OF_SCALARS_H;
						upper_index = i + j*NO_OF_SCALARS_H;
						upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
						lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
						total_volume = upper_volume + lower_volume;
						upper_weight = upper_volume/total_volume;
						lower_weight = lower_volume/total_volume;
						new_density_at_interface = upper_weight*state_new -> mass_densities[k*NO_OF_SCALARS + upper_index]
						+ lower_weight*state_new -> mass_densities[k*NO_OF_SCALARS + lower_index];
						vertical_flux_vector[j] = new_density_at_interface*vertical_flux_vector[j];
					}
				}
				// filling up the original vectors
				for (j = 0; j < NO_OF_LAYERS - 1; ++j)
				{
					a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
					c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
				}
				for (j = 0; j < NO_OF_LAYERS; ++j)
				{
					if (quantity_id == 1)
					{
						if (j == 0)
						{
							b_vector[j] = state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] - delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[0];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							b_vector[j] = state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j - 1];
						}
						else
						{
							b_vector[j] = state_new -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
						}
					}
					else
					{
						if (j == 0)
						{
							b_vector[j] = 1 - delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[0];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j - 1];
						}
						else
						{
							b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
						}
					}
					// the explicit component
					if (quantity_id == 0)
					{
						d_vector[j] =
						state_old -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
						+ delta_t*state_tendency -> mass_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
					}
					if (quantity_id == 1)
					{
						d_vector[j] =
						state_old -> entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
						+ delta_t*state_tendency -> entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
					}
					if (quantity_id == 2)
					{
						d_vector[j] =
						state_old -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
						+ delta_t*state_tendency -> condensed_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
					}
				}
				thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
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
						state_new -> entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] =
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

int thomas_algorithm(double a_vector[], double b_vector[], double c_vector[], double d_vector[], double c_prime_vector[], double d_prime_vector[], double solution_vector[], int solution_length)
{
	// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
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

int lu_5band_solver(double a_vector[], double b_vector[], double c_vector[], double d_vector[], double c_prime_vector[], double d_prime_vector[], double solution_vector[], int solution_length)
{
	for (int j = solution_length - 1; j >= 0; --j)
	{
		solution_vector[j] = 0;
	}
	return 0;
}

int sign(double x)
{
	if (x > 0)
	{
		return 1;
	}
	if (x < 0)
	{
		return -1;
	}
	return 0;
}





