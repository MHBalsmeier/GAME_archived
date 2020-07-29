/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "../../diagnostics/diagnostics.h"
#include <omp.h>
#include "../../spatial_operators/spatial_operators.h"
#include "atmostracers.h"

int thomas_algorithm(double [], double [], double [], double [], double [], double [], double [], int);
int sign(double);

int three_band_solver_hor_vel_adv(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Semi-implicit vertical advection of vertical momentum (Crank-Nicolson).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_velocity = malloc(NO_OF_LAYERS*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double delta_z;
	int i;
	for (i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		// diagnozing the vertical velocity
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			recov_hor_ver_pri(state_old -> velocity_gas, j, i, &vertical_velocity[j], grid);
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			if (j == NO_OF_LAYERS - 2)
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			else
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[j] = 0.5*vertical_velocity[j + 1]*delta_t/(2*delta_z);
			if (j == 0)
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			else
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[j] = -0.5*vertical_velocity[j]*delta_t/(2*delta_z);
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 + 0.5*delta_t*vertical_velocity[j]/(2*delta_z);
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				delta_z = grid -> z_vector[NO_OF_SCALARS_H + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 - 0.5*delta_t*vertical_velocity[j]/(2*delta_z);
			}
			else
			{
				b_vector[j] = 1;
			}
			d_vector[j] = state_old -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> velocity_gas[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_velocity);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_vel_adv(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Semi-implicit vertical advection of vertical momentum (Crank-Nicolson).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 2)*sizeof(double));
	double *b_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 2)*sizeof(double));
	double *d_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 2)*sizeof(double));
	double *d_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *vertical_velocity = malloc(NO_OF_LAYERS*sizeof(double));
	double *solution_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double delta_z;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// extracting the vertical velocity
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_velocity[j] = state_old -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
		}
		vertical_velocity[NO_OF_LAYERS - 1] = state_old -> velocity_gas[i + NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER];
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			delta_z = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 3)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[j] = 0.5*vertical_velocity[j + 1]*delta_t/delta_z;
		}
		delta_z = grid -> z_vector[i] - grid -> z_vector[2*NO_OF_VECTORS_PER_LAYER + i];
		c_vector[0] = -0.5*vertical_velocity[0]*delta_t/delta_z;
		for (int j = 1; j < NO_OF_LAYERS - 2; ++j)
		{
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[j] = -0.5*vertical_velocity[j]*delta_t/delta_z;
		}
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			b_vector[j] = 1;
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			if (j == 0)
			{
				d_vector[j] = state_old -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - 0.5*delta_t*vertical_velocity[j]*(-vertical_velocity[j + 1])/delta_z;
			}
			else if (j == NO_OF_LAYERS - 2)
			{
				d_vector[j] = state_old -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - 0.5*delta_t*vertical_velocity[j]*(vertical_velocity[j - 1] - 2*vertical_velocity[j + 1])/delta_z;
			}
			else
			{
				d_vector[j] = state_old -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - 0.5*delta_t*vertical_velocity[j]*(vertical_velocity[j - 1] - vertical_velocity[j + 1])/delta_z;
			}
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS - 1);
		state_new -> velocity_gas[i] = 0;
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_velocity);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_sound_waves(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, double delta_t, Grid *grid)
{
	double *a_vector = malloc((2*NO_OF_LAYERS - 2)*sizeof(double));
	double *b_vector = malloc((2*NO_OF_LAYERS - 1)*sizeof(double));
	double *c_vector = malloc((2*NO_OF_LAYERS - 2)*sizeof(double));
	double *d_vector = malloc((2*NO_OF_LAYERS - 1)*sizeof(double));
	double *c_prime_vector = malloc((2*NO_OF_LAYERS - 2)*sizeof(double));
	double *d_prime_vector = malloc((2*NO_OF_LAYERS - 1)*sizeof(double));
	double *solution_vector = malloc((2*NO_OF_LAYERS - 1)*sizeof(double));
	double *vertical_velocity = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *vertical_velocity_divergence = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_entropy_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_entropy_gradient = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *temp_interface_values = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *second_derivative_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double delta_z;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			vertical_entropy_vector[j] = diagnostics -> specific_entropy[i + j*NO_OF_SCALARS_H];
		}
		grad_v_scalar_column(vertical_entropy_vector, vertical_entropy_gradient, i, grid);
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_velocity[j] = state_new -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
		}
		for (int j = 1; j < NO_OF_LAYERS - 1; ++j)
		{
			second_derivative_vector[j] = state_old -> temp_gas[i + (j - 1)*NO_OF_SCALARS_H] - 2*state_old -> temp_gas[i + j*NO_OF_SCALARS_H] + state_old -> temp_gas[i + (j + 1)*NO_OF_SCALARS_H];
		}
		second_derivative_vector[0] = second_derivative_vector[1];
		second_derivative_vector[NO_OF_LAYERS - 1] = second_derivative_vector[NO_OF_LAYERS - 2];
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			temp_interface_values[j] = 0.5*(state_old -> temp_gas[i + j*NO_OF_SCALARS_H] + state_old -> temp_gas[i + (j + 1)*NO_OF_SCALARS_H]);
			// fourth order flux for stability
			temp_interface_values[j] += -1.0/12*(second_derivative_vector[j] + second_derivative_vector[j + 1]) - sign(state_old -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER])*1.0/12*(second_derivative_vector[j + 1] - second_derivative_vector[j]);
		}
		divv_v_columns(vertical_velocity, vertical_velocity_divergence, i, grid);
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			a_vector[2*j] = delta_t*C_D_V/C_D_P*C_D_P/delta_z - delta_t*C_D_V/C_D_P*0.5*vertical_entropy_gradient[j];
			delta_z = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[2*j + 1] = delta_t*((R_D/C_D_V - 1)*state_old -> temp_gas[i + (j + 1)*NO_OF_SCALARS_H] + temp_interface_values[j])/delta_z;
		}
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[2*j] = -delta_t*((R_D/C_D_V - 1)*state_old -> temp_gas[i + j*NO_OF_SCALARS_H] + temp_interface_values[j])/delta_z;
			delta_z = grid -> z_scalar[j*NO_OF_SCALARS_H + i] - grid -> z_scalar[(j + 1)*NO_OF_SCALARS_H + i];
			c_vector[2*j + 1] = -delta_t*C_D_V/C_D_P*C_D_P/delta_z - delta_t*C_D_V/C_D_P*0.5*vertical_entropy_gradient[j];
		}
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			b_vector[2*j] = 1 + delta_t*(R_D/C_D_V - 1)*diagnostics -> wind_field_divv_h[i + j*NO_OF_SCALARS_H];
			b_vector[2*j + 1] = 1;
			d_vector[2*j] = state_old -> temp_gas[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> temp_gas[j*NO_OF_SCALARS_H + i];
			d_vector[2*j + 1] = state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + delta_t*(-grid -> gravity_m[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + R_D/C_D_P*diagnostics -> pressure_gradient_1[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - R_D/C_D_P*diagnostics -> pressure_gradient_0_m[(j + 1)*NO_OF_VECTORS_PER_LAYER + i]);
		}
		b_vector[2*NO_OF_LAYERS - 2] =  1 + delta_t*(R_D/C_D_V - 1)*diagnostics -> wind_field_divv_h[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H];
		d_vector[2*NO_OF_LAYERS - 2] = state_old -> temp_gas[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> temp_gas[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, 2*NO_OF_LAYERS - 1);
		// writing the result into the new state
		state_new -> velocity_gas[i] = 0;
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			state_new -> temp_gas[i + j*NO_OF_SCALARS_H] = solution_vector[2*j];
			state_new -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[2*j + 1];
		}
		state_new -> temp_gas[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H] = solution_vector[2*(NO_OF_LAYERS - 1)];
	}
	free(vertical_entropy_vector);
	free(vertical_entropy_gradient);
	free(second_derivative_vector);
	free(temp_interface_values);
	free(vertical_velocity);
	free(vertical_velocity_divergence);
	free(solution_vector);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_den_dry(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of dry mass (Euler).
	Procedure derived in Kompendium.
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_flux_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	int i;
	double area;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// diagnozing the vertical flux
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				vertical_contravariant_normalized(state_new -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
			else
				vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
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
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> density_dry[j*NO_OF_SCALARS_H + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_flux_vector);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_entropy_density_gas(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of the entropy of the gas phase (Euler).
	Procedure derived in Kompendium.
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_flux_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double area;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// diagnozing the vertical flux
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				vertical_contravariant_normalized(state_new -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
			else
				vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
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
			d_vector[j] = state_old -> entropy_density_gas[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> entropy_density_gas[j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> entropy_density_gas[j*NO_OF_SCALARS_H + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_flux_vector);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_tracers(State *state_old, State *state_new, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of tracers (Euler).
	Procedure derived in Kompendium.
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_flux_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double area;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		for (int k = 0; k < NO_OF_CONDENSATED_TRACERS; ++k)
		{
			// tracer masses
			// diagnozing the vertical flux
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				{
					vertical_contravariant_normalized(state_new -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
					// small metric terms neglected here
					vertical_flux_vector[j] -= ret_sink_velocity(i, 0, 0.001);
				}
				else
					vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] = area*vertical_flux_vector[j];
			}
			// filling up the original vectors
			for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
				c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			}
			for (int j = 0; j < NO_OF_LAYERS; ++j)
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
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_new -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_new -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
					state_new -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
			}
			// tracer temberature densities
			// diagnozing the vertical flux
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				{
					vertical_contravariant_normalized(state_new -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
					// small metric term negelcted here
					vertical_flux_vector[j] -= ret_sink_velocity(i, 0, 0.001);
				}
				else
					vertical_flux_vector[j] = state_new -> velocity_gas[k*NO_OF_SCALARS + NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				vertical_flux_vector[j] = area*vertical_flux_vector[j];
			}
			// filling up the original vectors
			for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
				c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			}
			for (int j = 0; j < NO_OF_LAYERS; ++j)
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
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_new -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_new -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] < 0)
					state_new -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] = 0;
			}
		}
		// water vaour
		// diagnozing the vertical flux
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				vertical_contravariant_normalized(state_new -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
			else
				vertical_flux_vector[j] = state_new -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
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
			d_vector[j] = state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_new -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
			if (state_new -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
				state_new -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
		}
	}
	free(solution_vector);
	free(vertical_flux_vector);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int thomas_algorithm(double a_vector[], double b_vector[], double c_vector[], double d_vector[], double c_prime_vector[], double d_prime_vector[], double solution_vector[], int solution_length)
{
	// https://de.wikipedia.org/wiki/Thomas-Algorithmus
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




