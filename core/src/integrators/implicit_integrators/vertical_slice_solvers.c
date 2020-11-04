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
int sign(double);

int three_band_solver_hor_vel_adv(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Semi-implicit vertical advection of vertical momentum (Crank-Nicolson).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double delta_z;
	int i, j;
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
			recov_hor_ver_pri(state_old -> velocity_gas, j, i, &vertical_velocity[j], grid);
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
				- 0.5*delta_t*vertical_velocity[j]
				*(state_old -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - state_old -> velocity_gas[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i])
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
				- 0.5*delta_t*vertical_velocity[j]*state_old -> velocity_gas[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i]/delta_z;
			}
			else
			{
				b_vector[j] = 1;
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i]- grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
				// right hand side
				// old value
				d_vector[j] = state_old -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i]
				// the old time step part of Crank-Nicolson
				- 0.5*delta_t*vertical_velocity[j]
				*(state_old -> velocity_gas[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - state_old -> velocity_gas[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i])/delta_z;
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
			a_vector[j] = 0.5*vertical_velocity[j + 1]*delta_t/(2*delta_z);
			if (j == 0)
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			}
			else
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			}
			c_vector[j] = -0.5*vertical_velocity[j]*delta_t/(2*delta_z);
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

int three_band_solver_ver_vel_adv(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Semi-implicit vertical advection of vertical momentum (Crank-Nicolson).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double delta_z;
	int i, j;
	#pragma omp parallel for private(delta_z, j)
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{	
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[NO_OF_LAYERS - 2];
		double b_vector[NO_OF_LAYERS - 1];
		double c_vector[NO_OF_LAYERS - 2];
		double d_vector[NO_OF_LAYERS - 1];
		double c_prime_vector[NO_OF_LAYERS - 2];
		double d_prime_vector[NO_OF_LAYERS - 1];
		double vertical_velocity[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS - 1];
		// extracting the vertical velocity
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_velocity[j] = state_old -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
		}
		// filling up b vector and d vector
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			if (j == 0)
			{
				delta_z = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 + delta_t*vertical_velocity[j]/(2*delta_z);
				// right hand side
				// old value
				d_vector[j] = vertical_velocity[j]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]
				// the old time step part of Crank-Nicolson
				- 0.5*delta_t*vertical_velocity[j]*(vertical_velocity[j] - vertical_velocity[j + 1])/delta_z;
			}
			else if (j == NO_OF_LAYERS - 2)
			{
				b_vector[j] = 1;
				delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
				// right hand side
				// old value
				d_vector[j] = vertical_velocity[j]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]
				// the old time step part of Crank-Nicolson
				- 0.5*delta_t*vertical_velocity[j]*vertical_velocity[j - 1]/delta_z;
			}
			else
			{
				b_vector[j] = 1;
				delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[i + (j + 2)*NO_OF_VECTORS_PER_LAYER];
				// right hand side
				// old value
				d_vector[j] = vertical_velocity[j]
				// explicit part
				+ delta_t*state_tendency -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]
				// the old time step part of Crank-Nicolson
				- 0.5*delta_t*vertical_velocity[j]*(vertical_velocity[j - 1] - vertical_velocity[j + 1])/delta_z;
			}
		}
		// filling up the a vector
		for (j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			delta_z = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 3)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[j] = 0.5*vertical_velocity[j + 1]*delta_t/delta_z;
		}
		// filling up the c vector
		delta_z = grid -> z_vector[NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[2*NO_OF_VECTORS_PER_LAYER + i];
		c_vector[0] = -0.5*vertical_velocity[0]*delta_t/delta_z;
		for (j = 1; j < NO_OF_LAYERS - 2; ++j)
		{
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[j] = -0.5*vertical_velocity[j]*delta_t/delta_z;
		}
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS - 1);
		// writing the solution into the new state
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			state_new -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] = solution_vector[j];
		}
	}
	return 0;
}

int three_band_solver_ver_sound_waves(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	double delta_z, upper_weight, lower_weight, upper_volume, lower_volume, total_volume, damping_coeff, damping_coeff_max, damping_start_height, z_above_damping, damping_start_height_over_toa, c_h_v, c_h_p, R_h, density_d_value, density_v_value;
	// This is for Klemp (2008).
	get_damping_layer_properties(&damping_start_height_over_toa, &damping_coeff_max);
	damping_start_height = damping_start_height_over_toa*grid -> z_vector[0];
	int upper_index, lower_index;
	int j;
	#pragma omp parallel for private(upper_index, lower_index, delta_z, upper_weight, lower_weight, upper_volume, lower_volume, total_volume, damping_coeff, z_above_damping, j, c_h_v, c_h_p, R_h, density_d_value, density_v_value)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[2*NO_OF_LAYERS - 2];
		double b_vector[2*NO_OF_LAYERS - 1];
		double c_vector[2*NO_OF_LAYERS - 2];
		double d_vector[2*NO_OF_LAYERS - 1];
		double c_prime_vector[2*NO_OF_LAYERS - 2];
		double d_prime_vector[2*NO_OF_LAYERS - 1]; 
		double solution_vector[2*NO_OF_LAYERS - 1];
		double vertical_entropy_vector[NO_OF_LAYERS];
		double vertical_entropy_gradient[NO_OF_LAYERS - 1];
		double temp_interface_values[NO_OF_LAYERS - 1];
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			vertical_entropy_vector[j] = diagnostics -> specific_entropy_dry[i + j*NO_OF_SCALARS_H];
		}
		grad_v_scalar_column(vertical_entropy_vector, vertical_entropy_gradient, i, grid);
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
            lower_index = i + (j + 1)*NO_OF_SCALARS_H;
            upper_index = i + j*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            upper_weight = upper_volume/total_volume;
            lower_weight = lower_volume/total_volume;
			temp_interface_values[j] = upper_weight*state_old -> temperature_gas[i + j*NO_OF_SCALARS_H] + lower_weight*state_old -> temperature_gas[i + (j + 1)*NO_OF_SCALARS_H];
		}
		// filling up the original vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			density_d_value = state_old -> density_dry[j*NO_OF_SCALARS_H + i];
			density_v_value = state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			c_h_v = spec_heat_cap_diagnostics_v(density_d_value, density_v_value);
			c_h_p = spec_heat_cap_diagnostics_p(density_d_value, density_v_value);
			R_h = gas_constant_diagnostics(density_d_value, density_v_value);
			
			// determining the elements of a_vector
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			a_vector[2*j] = delta_t*c_h_v/c_h_p*c_h_p/delta_z - delta_t*c_h_v/c_h_p*0.5*vertical_entropy_gradient[j];
			delta_z = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[2*j + 1] = delta_t*((R_h/c_h_v - 1)*state_old -> temperature_gas[i + (j + 1)*NO_OF_SCALARS_H] + temp_interface_values[j])/delta_z;
			
			// determining the elements of c_vector
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[2*j] = -delta_t*((R_h/c_h_v - 1)*state_old -> temperature_gas[i + j*NO_OF_SCALARS_H] + temp_interface_values[j])/delta_z;
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			c_vector[2*j + 1] = -delta_t*c_h_v/c_h_p*c_h_p/delta_z - delta_t*c_h_v/c_h_p*0.5*vertical_entropy_gradient[j];
		}
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			b_vector[2*j] = 1;
			b_vector[2*j + 1] = 1;
			d_vector[2*j] = diagnostics -> temperature_gas_explicit[j*NO_OF_SCALARS_H + i];
			d_vector[2*j + 1] = state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
		}
		b_vector[2*NO_OF_LAYERS - 2] =  1;
		d_vector[2*NO_OF_LAYERS - 2] = diagnostics -> temperature_gas_explicit[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, 2*NO_OF_LAYERS - 1);
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
			state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[2*j + 1]/(1 + delta_t*damping_coeff);
		}
	}
	return 0;
}

int three_band_solver_ver_den_dry(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of dry mass (Euler).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double area;
	int j;
	#pragma omp parallel for private(area, j)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[NO_OF_LAYERS - 1];
		double b_vector[NO_OF_LAYERS];
		double c_vector[NO_OF_LAYERS - 1];
		double d_vector[NO_OF_LAYERS];
		double c_prime_vector[NO_OF_LAYERS - 1];
		double d_prime_vector[NO_OF_LAYERS];
		double vertical_flux_vector[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS];
		// diagnozing the vertical flux
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (j = 0; j < NO_OF_LAYERS; ++j)
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
			d_vector[j] = state_old -> density_dry[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> density_dry[j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> density_dry[j*NO_OF_SCALARS_H + i] = solution_vector[j];
		}
	}
	return 0;
}

int three_band_solver_ver_entropy_density_dry(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of the entropy of the dry gas phase (Euler).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double area;
	int j;
	#pragma omp parallel for private(area, j)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[NO_OF_LAYERS - 1];
		double b_vector[NO_OF_LAYERS];
		double c_vector[NO_OF_LAYERS - 1];
		double d_vector[NO_OF_LAYERS];
		double c_prime_vector[NO_OF_LAYERS - 1];
		double d_prime_vector[NO_OF_LAYERS];
		double vertical_flux_vector[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS];
		// diagnozing the vertical flux
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (j = 0; j < NO_OF_LAYERS; ++j)
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
			d_vector[j] = state_old -> entropy_density_dry[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> entropy_density_dry[j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> entropy_density_dry[j*NO_OF_SCALARS_H + i] = solution_vector[j];
		}
	}
	return 0;
}

int three_band_solver_ver_entropy_density_gaseous_tracers(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of the entropy of the gaseous tracers (Euler).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double area;
	int j, k;
	#pragma omp parallel for private(area, j, k)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[NO_OF_LAYERS - 1];
		double b_vector[NO_OF_LAYERS];
		double c_vector[NO_OF_LAYERS - 1];
		double d_vector[NO_OF_LAYERS];
		double c_prime_vector[NO_OF_LAYERS - 1];
		double d_prime_vector[NO_OF_LAYERS];
		double vertical_flux_vector[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS];
		for (k = 0; k < NO_OF_GASEOUS_TRACERS; ++k)
		{
			// diagnozing the vertical flux
			for (j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] = area*vertical_flux_vector[j];
			}
			// filling up the original vectors
			for (j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
				c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			}
			for (j = 0; j < NO_OF_LAYERS; ++j)
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
				d_vector[j] = state_old -> tracer_entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i]
				+ delta_t*state_tendency -> tracer_entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			}
			thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_new -> tracer_entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_new -> tracer_entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
				{
					state_new -> tracer_entropy_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
				}
			}
		}
	}
	return 0;
}

int three_band_solver_ver_tracers(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of tracers (Euler).
	procedure derived in https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf
	*/
	double area;
	int j, k;
	#pragma omp parallel for private(area, j, k)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// for meanings of these vectors look into the definition of the function thomas_algorithm
		double a_vector[NO_OF_LAYERS - 1];
		double b_vector[NO_OF_LAYERS];
		double c_vector[NO_OF_LAYERS - 1];
		double d_vector[NO_OF_LAYERS];
		double c_prime_vector[NO_OF_LAYERS - 1];
		double d_prime_vector[NO_OF_LAYERS];
		double solution_vector[NO_OF_LAYERS];
		double vertical_flux_vector[NO_OF_LAYERS];
		for (k = 0; k < NO_OF_CONDENSED_TRACERS; ++k)
		{
			// determining the vertical flux
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] -= ret_sink_velocity(i, 0, 0.001);
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] = area*vertical_flux_vector[j];
			}
			// filling up the original vectors
			for (j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
				c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			}
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				if (j == 0)
				{
					b_vector[j] = 1 - delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[0];
				}
				else
				{
					b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
				}
				d_vector[j] = state_old -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			}
			thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_new -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_new -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
				{
					state_new -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
				}
			}
			// tracer temperature densities
			// diagnozing the vertical flux
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				vertical_flux_vector[j] = state_new -> velocity_gas[k*NO_OF_SCALARS + NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] -= ret_sink_velocity(i, 0, 0.001);
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] = area*vertical_flux_vector[j];
			}
			// filling up the original vectors
			for (j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
				c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			}
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				if (j == 0)
				{
					b_vector[j] = 1 - delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[0];
				}
				else
				{
					b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
				}
				d_vector[j] = state_old -> tracer_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			}
			thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
			for (j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_new -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_new -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] < 0)
				{
					state_new -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] = 0;
				}
			}
		}
		// water vaour
		// diagnozing the vertical flux
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (j = 0; j < NO_OF_LAYERS; ++j)
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
			d_vector[j] = state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
			// limiter
			if (state_new -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
			{
				state_new -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
			}
		}
	}
	return 0;
}

int thomas_algorithm(double a_vector[], double b_vector[], double c_vector[], double d_vector[], double c_prime_vector[], double d_prime_vector[], double solution_vector[], int solution_length)
{
	// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	c_prime_vector[0] = c_vector[0]/b_vector[0];
	for (int j = 1; j < solution_length - 1; ++j)
	{
		c_prime_vector[j] = c_vector[j]/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1]);
	}
	d_prime_vector[0] = d_vector[0]/b_vector[0];
	for (int j = 1; j < solution_length; ++j)
	{
		d_prime_vector[j] = (d_vector[j] - d_prime_vector[j - 1]*a_vector[j - 1])/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1]);
	}
	solution_vector[solution_length - 1] = d_prime_vector[solution_length - 1];
	for (int j = solution_length - 2; j >= 0; --j)
	{
		solution_vector[j] = d_prime_vector[j] - c_prime_vector[j]*solution_vector[j + 1];
	}
	return 0;
}

int sign(double x)
{
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}





