/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains the implicit vertical solvers.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"
#include "../constituents/constituents.h"
#include "../subgrid_scale/subgrid_scale.h"

int thomas_algorithm(double [], double [], double [], double [], double [], int);

int three_band_solver_ver_waves(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings,
Config *config, double delta_t, Grid *grid, int rk_step)
{
	/*
	This is the implicit vertical solver for the main fluid constituent.
	*/
	
	// declaring and defining some variables that will be needed later on
	int lower_index, base_index;
	double impl_weight = config -> impl_thermo_weight;
	double c_v = spec_heat_capacities_v_gas(0);
	double c_p = spec_heat_capacities_p_gas(0);
	double r_d = specific_gas_constants(0);
	// This is for Klemp (2008).
	double damping_coeff, damping_start_height, z_above_damping, temperature_gas_lowest_layer_old, temperature_gas_lowest_layer_new;
	damping_start_height = config -> damping_start_height_over_toa*grid -> z_vector[0];
	
	// partial derivatives new time step weight
	double partial_deriv_new_time_step_weight = 0.5;
	
	int gas_phase_first_index = NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS;
	
	// calculating the sensible power flux density
	if (config -> soil_on == 1)
	{
		#pragma omp parallel for private(base_index, temperature_gas_lowest_layer_old, temperature_gas_lowest_layer_new)
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			base_index = NO_OF_SCALARS - NO_OF_SCALARS_H + i;
			
			// gas temperature in the lowest layer
			temperature_gas_lowest_layer_old = (grid -> exner_bg[base_index] + state_old -> exner_pert[base_index])
			*(grid -> theta_bg[base_index] + state_old -> theta_pert[base_index]);
			temperature_gas_lowest_layer_new = (grid -> exner_bg[base_index] + state_new -> exner_pert[base_index])
			*(grid -> theta_bg[base_index] + state_new -> theta_pert[base_index]);
			
			// the sensible power flux density
			diagnostics -> power_flux_density_sensible[i] = 0.5*c_v*(state_new -> rho[gas_phase_first_index + base_index]
			*(temperature_gas_lowest_layer_old - state_old -> temperature_soil[i])
			+ state_old -> rho[gas_phase_first_index + base_index]
			*(temperature_gas_lowest_layer_new - state_new -> temperature_soil[i]))/diagnostics -> scalar_flux_resistance[i];
			
			// contribution of sensible heat to rhotheta
			state_tendency -> rhotheta[base_index] 
			+= -grid -> area[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i]*diagnostics -> power_flux_density_sensible[i]
			/((grid -> exner_bg[base_index] + state_new -> exner_pert[base_index])*c_p)/grid -> volume[base_index];
		}
	}
	
	// loop over all columns
	#pragma omp parallel for private(lower_index, damping_coeff, z_above_damping, base_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		
		// for meanings of these vectors look into the Kompendium
		double c_vector[NO_OF_LAYERS - 2 + config -> soil_on*NO_OF_SOIL_LAYERS];
		double d_vector[NO_OF_LAYERS - 1 + config -> soil_on*NO_OF_SOIL_LAYERS];
		double e_vector[NO_OF_LAYERS - 2 + config -> soil_on*NO_OF_SOIL_LAYERS];
		double r_vector[NO_OF_LAYERS - 1 + config -> soil_on*NO_OF_SOIL_LAYERS];
		double rho_expl[NO_OF_LAYERS];
		double rhotheta_expl[NO_OF_LAYERS];
		double theta_pert_expl[NO_OF_LAYERS];
		double exner_pert_expl[NO_OF_LAYERS];
		double theta_int_new[NO_OF_LAYERS - 1];
		double solution_vector[NO_OF_LAYERS - 1 + config -> soil_on*NO_OF_SOIL_LAYERS];
		double rho_int_old[NO_OF_LAYERS - 1];
		double rho_int_expl[NO_OF_LAYERS - 1];
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
			base_index = i + j*NO_OF_SCALARS_H;
			// explicit density
			rho_expl[j] = state_old -> rho[gas_phase_first_index + base_index]
			+ delta_t*state_tendency -> rho[gas_phase_first_index + base_index];
			// explicit potential temperature density
			rhotheta_expl[j] = state_old -> rhotheta[base_index] + delta_t*state_tendency -> rhotheta[base_index];
			if (rk_step == 0)
			{
				// old time step partial derivatives of theta and Pi (divided by the volume)
				alpha[j] = -state_old -> rhotheta[base_index]/pow(state_old -> rho[gas_phase_first_index + base_index], 2)
				/grid -> volume[base_index];
				beta[j] = 1.0/state_old -> rho[gas_phase_first_index + base_index]/grid -> volume[base_index];
				gamma[j] = r_d/(c_v*state_old -> rhotheta[base_index])
				*(grid -> exner_bg[base_index] + state_old -> exner_pert[base_index])/grid -> volume[base_index];
			}
			else
			{
				// old time step partial derivatives of theta and Pi
				alpha_old[j] = -state_old -> rhotheta[base_index]/pow(state_old -> rho[gas_phase_first_index + base_index], 2);
				beta_old[j] = 1.0/state_old -> rho[gas_phase_first_index + base_index];
				gamma_old[j] = r_d/(c_v*state_old -> rhotheta[base_index])*(grid -> exner_bg[base_index] + state_old -> exner_pert[base_index]);
				// new time step partial derivatives of theta and Pi
				alpha_new[j] = -state_new -> rhotheta[base_index]/pow(state_new -> rho[gas_phase_first_index + base_index], 2);
				beta_new[j] = 1.0/state_new -> rho[gas_phase_first_index + base_index];
				gamma_new[j] = r_d/(c_v*state_new -> rhotheta[base_index])*(grid -> exner_bg[base_index] + state_new -> exner_pert[base_index]);
				// interpolation in time and dividing by the volume
				alpha[j] = ((1.0 - partial_deriv_new_time_step_weight)*alpha_old[j] + partial_deriv_new_time_step_weight*alpha_new[j])/grid -> volume[base_index];
				beta[j] = ((1.0 - partial_deriv_new_time_step_weight)*beta_old[j] + partial_deriv_new_time_step_weight*beta_new[j])/grid -> volume[base_index];
				gamma[j] = ((1.0 - partial_deriv_new_time_step_weight)*gamma_old[j] + partial_deriv_new_time_step_weight*gamma_new[j])/grid -> volume[base_index];
			}
			// explicit potential temperature perturbation
			theta_pert_expl[j] = state_old -> theta_pert[base_index] + delta_t*grid -> volume[base_index]*(
			alpha[j]*state_tendency -> rho[gas_phase_first_index + base_index] + beta[j]*state_tendency -> rhotheta[base_index]);
			// explicit Exner pressure perturbation
			exner_pert_expl[j] = state_old -> exner_pert[base_index] + delta_t*grid -> volume[base_index]*gamma[j]*state_tendency -> rhotheta[base_index];
		}
		
		// determining the interface values
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			lower_index = i + (j + 1)*NO_OF_SCALARS_H;
			rho_int_old[j] = 0.5*(state_old -> rho[gas_phase_first_index + base_index] + state_old -> rho[gas_phase_first_index + lower_index]);
			rho_int_expl[j] = 0.5*(rho_expl[j] + rho_expl[j + 1]);
			theta_int_new[j] = 0.5*(state_new -> rhotheta[base_index]/state_new -> rho[gas_phase_first_index + base_index]
			+ state_new -> rhotheta[lower_index]/state_new -> rho[gas_phase_first_index + lower_index]);
		}
		
		// filling up the coefficient vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			lower_index = i + (j + 1)*NO_OF_SCALARS_H;
			// main diagonal
			d_vector[j] = -pow(theta_int_new[j], 2)*(gamma[j] + gamma[j + 1])
			+ 0.5*(grid -> exner_bg[base_index] - grid -> exner_bg[lower_index])
			*(alpha[j + 1] - alpha[j] + theta_int_new[j]*(beta[j + 1] - beta[j]))
			- (grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])/(impl_weight*pow(delta_t, 2)*c_p*rho_int_old[j])
			*(2.0/grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] + delta_t*state_old -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]*0.5
			*(-1.0/grid -> volume[base_index] + 1.0/grid -> volume[lower_index]));
			// right hand side
			r_vector[j] = -(state_old -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] + delta_t*state_tendency -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER])
			*(grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])
			/(impl_weight*pow(delta_t, 2)*c_p)
			+ theta_int_new[j]*(exner_pert_expl[j] - exner_pert_expl[j + 1])/delta_t
			+ 0.5/delta_t*(theta_pert_expl[j] + theta_pert_expl[j + 1])*(grid -> exner_bg[base_index] - grid -> exner_bg[lower_index])
			- (grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])/(impl_weight*pow(delta_t, 2)*c_p)
			*state_old -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]*rho_int_expl[j]/rho_int_old[j];
		}
		for (int j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			lower_index = i + (j + 1)*NO_OF_SCALARS_H;
			// lower diagonal
			c_vector[j] = theta_int_new[j + 1]*gamma[j + 1]*theta_int_new[j]
			+ 0.5*(grid -> exner_bg[lower_index] - grid -> exner_bg[(j + 2)*NO_OF_SCALARS_H + i])
			*(alpha[j + 1] + beta[j + 1]*theta_int_new[j])
			- (grid -> z_scalar[lower_index] - grid -> z_scalar[(j + 2)*NO_OF_SCALARS_H + i])/(impl_weight*delta_t*c_p)*0.5
			*state_old -> wind[i + (j + 2)*NO_OF_VECTORS_PER_LAYER]/(grid -> volume[lower_index]*rho_int_old[j + 1]);
			// upper diagonal
			e_vector[j] = theta_int_new[j]*gamma[j + 1]*theta_int_new[j + 1]
			- 0.5*(grid -> exner_bg[base_index] - grid -> exner_bg[lower_index])
			*(alpha[j + 1] + beta[j + 1]*theta_int_new[j + 1])
			+ (grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])/(impl_weight*delta_t*c_p)*0.5
			*state_old -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]/(grid -> volume[lower_index]*rho_int_old[j]);
		}
		
		// soil components of the matrix
		if (config -> soil_on == 1)
		{
			// calculating the explicit part of the heat flux density
			double heat_flux_density_expl[NO_OF_SOIL_LAYERS];
			for (int j = 0; j < NO_OF_SOIL_LAYERS - 1; ++j)
			{
				heat_flux_density_expl[j]
				= -grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(state_old -> temperature_soil[i + j*NO_OF_SCALARS_H]
				- state_old -> temperature_soil[i + (j + 1)*NO_OF_SCALARS_H])
				/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
			}
			heat_flux_density_expl[NO_OF_SOIL_LAYERS - 1]
			= -grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(state_old -> temperature_soil[i + (NO_OF_SOIL_LAYERS - 1)*NO_OF_SCALARS_H]
			- grid -> t_const_soil[i])
			/(2*(grid -> z_soil_center[NO_OF_SOIL_LAYERS - 1] - grid -> z_t_const));
			
			// calculating the explicit part of the temperature change
			r_vector[NO_OF_LAYERS - 1]
			// old temperature
			= state_old -> temperature_soil[i]
			// sensible heat flux
			+ (diagnostics -> power_flux_density_sensible[i]
			// latent heat flux
			+ diagnostics -> power_flux_density_latent[i]
			// shortwave inbound radiation
			+ forcings -> sfc_sw_in[i]
			// longwave outbound radiation
			- forcings -> sfc_lw_out[i]
			// heat conduction from below
			+ 0.5*heat_flux_density_expl[0])
			/((grid -> z_soil_interface[0] - grid -> z_soil_interface[1])*grid -> sfc_rho_c[i])*delta_t;
			
			// loop over all soil layers below the first layer
			for (int j = 1; j < NO_OF_SOIL_LAYERS; ++j)
			{
				
				r_vector[j + NO_OF_LAYERS - 1]
				// old temperature
				= state_old -> temperature_soil[i + j*NO_OF_SCALARS_H]
				// heat conduction from above
				+ 0.5*(-heat_flux_density_expl[j - 1]
				// heat conduction from below
				+ heat_flux_density_expl[j])
				/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])*delta_t;
			}
			
			// the diagonal component
			for (int j = 0; j < NO_OF_SOIL_LAYERS; ++j)
			{
				if (j == 0)
				{
					d_vector[j + NO_OF_LAYERS - 1] = 1.0 + 0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
					/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
					*1.0/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
				}
				else if (j == NO_OF_SOIL_LAYERS - 1)
				{
					d_vector[j + NO_OF_LAYERS - 1] = 1.0 + 0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
					/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
					*1.0/(grid -> z_soil_center[j - 1] - grid -> z_soil_center[j]);
				}
				else
				{
					d_vector[j + NO_OF_LAYERS - 1] = 1.0 + 0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
					/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
					*(1.0/(grid -> z_soil_center[j - 1] - grid -> z_soil_center[j])
					+ 1.0/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]));
				}
			}
			// the off-diagonal components
			c_vector[NO_OF_LAYERS - 2] = 0.0;
			e_vector[NO_OF_LAYERS - 2] = 0.0;
			for (int j = 0; j < NO_OF_SOIL_LAYERS - 1; ++j)
			{
				c_vector[j + NO_OF_LAYERS - 1] = -0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
				/((grid -> z_soil_interface[j + 1] - grid -> z_soil_interface[j + 2])*grid -> sfc_rho_c[i])
				/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
				e_vector[j + NO_OF_LAYERS - 1] = -0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
				/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
				/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
			}
		}
		
		
		// calling the algorithm to solve the system of linear equations
		thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_LAYERS - 1 + config -> soil_on*NO_OF_SOIL_LAYERS);
		
		// Klemp (2008) upper boundary layer
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			z_above_damping = grid -> z_vector[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] - damping_start_height;
			if (z_above_damping < 0.0)
			{
				damping_coeff = 0.0;
			}
			else
			{
				damping_coeff = config -> damping_coeff_max*pow(sin(0.5*M_PI*z_above_damping/(grid -> z_vector[0] - damping_start_height)), 2);
			}
			solution_vector[j] = solution_vector[j]/(1.0 + delta_t*damping_coeff);
		}
		
		/*
		Writing the result into the new state.
		--------------------------------------
		*/
		// mass density
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			if (j == 0)
			{
				state_new -> rho[gas_phase_first_index + base_index]
				= rho_expl[j] + delta_t*(solution_vector[j])/grid -> volume[base_index];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				state_new -> rho[gas_phase_first_index + base_index]
				= rho_expl[j] + delta_t*(-solution_vector[j - 1])/grid -> volume[base_index];
			}
			else
			{
				state_new -> rho[gas_phase_first_index + base_index]
				= rho_expl[j] + delta_t*(-solution_vector[j - 1] + solution_vector[j])/grid -> volume[base_index];
			}
		}
		// potential temperature density
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			if (j == 0)
			{
				state_new -> rhotheta[base_index]
				= rhotheta_expl[j] + delta_t*(theta_int_new[j]*solution_vector[j])/grid -> volume[base_index];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				state_new -> rhotheta[base_index]
				= rhotheta_expl[j] + delta_t*(-theta_int_new[j - 1]*solution_vector[j - 1])/grid -> volume[base_index];
			}
			else
			{
				state_new -> rhotheta[base_index]
				= rhotheta_expl[j] + delta_t*(-theta_int_new[j - 1]*solution_vector[j - 1] + theta_int_new[j]*solution_vector[j])
				/grid -> volume[base_index];
			}
		}
		// vertical velocity
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			density_interface_new
			= 0.5*(state_new -> rho[gas_phase_first_index + base_index]
			+ state_new -> rho[gas_phase_first_index + i + (j + 1)*NO_OF_SCALARS_H]);
			state_new -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]
			= (2.0*solution_vector[j]/grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] - density_interface_new*state_old -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER])
			/rho_int_old[j];
		}
		// potential temperature perturbation
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			state_new -> theta_pert[base_index] = state_new -> rhotheta[base_index]
			/state_new -> rho[gas_phase_first_index + base_index]
			- grid -> theta_bg[base_index];
		}
		// Exner pressure perturbation
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			base_index = i + j*NO_OF_SCALARS_H;
			state_new -> exner_pert[base_index] = state_old -> exner_pert[base_index] + grid -> volume[base_index]
			*gamma[j]*(state_new -> rhotheta[base_index] - state_old -> rhotheta[base_index]);
		}
		
		// soil temperature
		if (config -> soil_on == 1)
		{
			for (int j = 0; j < NO_OF_SOIL_LAYERS; ++j)
			{
				state_new -> temperature_soil[i + j*NO_OF_SCALARS_H] = solution_vector[NO_OF_LAYERS - 1 + j];
			}
		}
		
	} // end of the column (index i) loop
	return 0;
}

int three_band_solver_gen_densities(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config *config, double delta_t, Grid *grid)
{
	// Vertical advection of generalized densities (of tracers) with 3-band matrices.
	// mass densities, density x temperatures
	int no_of_relevant_constituents;
	double impl_weight, expl_weight;
	impl_weight = 0.5;
	expl_weight = 1.0 - impl_weight;
	for (int quantity_id = 0; quantity_id < 2; ++quantity_id)
	{
		no_of_relevant_constituents = 0;
		// mass densities
		if (quantity_id == 0)
		{
			// all constituents have a mass density
			no_of_relevant_constituents = NO_OF_CONSTITUENTS; // the main gaseous constituent is excluded later
		}
		// density x temperature fields
		if (quantity_id == 1)
		{
			// in this case, all the condensed constituents have a density x temperature field
			if(config -> assume_lte == 0)
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
					int lower_index, upper_index, base_index;
					
					// diagnozing the vertical fluxes
					for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
					{
						base_index = i + j*NO_OF_SCALARS_H;
						vertical_flux_vector_impl[j] = state_old -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
						vertical_flux_vector_rhs[j] = state_new -> wind[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
						// preparing the vertical interpolation
						lower_index = i + (j + 1)*NO_OF_SCALARS_H;
						upper_index = base_index;
						// For condensed constituents, a sink velocity must be added.
						// precipitation
						// snow
						if (k < NO_OF_CONDENSED_CONSTITUENTS/4)
						{
							vertical_flux_vector_impl[j] -= config -> snow_velocity;
							vertical_flux_vector_rhs[j] -= config -> snow_velocity;
						}
						// rain
						else if (k < NO_OF_CONDENSED_CONSTITUENTS/2)
						{
							vertical_flux_vector_impl[j] -= config -> rain_velocity;
							vertical_flux_vector_rhs[j] -= config -> rain_velocity;
						}
						// clouds
						else if (k < NO_OF_CONDENSED_CONSTITUENTS)
						{
							vertical_flux_vector_impl[j] -= config -> cloud_droplets_velocity;
							vertical_flux_vector_rhs[j] -= config -> cloud_droplets_velocity;
						}
						// multiplying the vertical velocity by the area
						area = grid -> area[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
						vertical_flux_vector_impl[j] = area*vertical_flux_vector_impl[j];
						vertical_flux_vector_rhs[j] = area*vertical_flux_vector_rhs[j];
						// old density at the interface
						if (vertical_flux_vector_rhs[j] >= 0.0)
						{
							density_old_at_interface = state_old -> rho[k*NO_OF_SCALARS + lower_index];
						}
						else
						{
							density_old_at_interface = state_old -> rho[k*NO_OF_SCALARS + upper_index];
						}
						vertical_flux_vector_rhs[j] = density_old_at_interface*vertical_flux_vector_rhs[j];
					}
					
					/*
					Now we proceed to solving the vertical tridiagonal problems.
					*/
					// filling up the original vectors
					for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
					{
						base_index = i + j*NO_OF_SCALARS_H;
						if (vertical_flux_vector_impl[j] >= 0.0)
						{
							c_vector[j] = 0.0;
							e_vector[j] = -impl_weight*delta_t/grid -> volume[base_index]*vertical_flux_vector_impl[j];
						}
						else
						{
							c_vector[j] = impl_weight*delta_t/grid -> volume[i + (j + 1)*NO_OF_SCALARS_H]*vertical_flux_vector_impl[j];
							e_vector[j] = 0.0;
						}
					}
					for (int j = 0; j < NO_OF_LAYERS; ++j)
					{
						base_index = i + j*NO_OF_SCALARS_H;
						if (j == 0)
						{
							if (vertical_flux_vector_impl[0] >= 0.0)
							{
								d_vector[j] = 1.0;
							}
							else
							{
								d_vector[j] = 1.0 - impl_weight*delta_t/grid -> volume[base_index]*vertical_flux_vector_impl[0];
							}
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							if (vertical_flux_vector_impl[j - 1] >= 0.0)
							{
								d_vector[j] = 1.0 + impl_weight*delta_t/grid -> volume[base_index]*vertical_flux_vector_impl[j - 1];
							}
							else
							{
								d_vector[j] = 1.0;
							}
							// precipitation
							// snow
							if (k < NO_OF_CONDENSED_CONSTITUENTS/4)
							{
								d_vector[j] += impl_weight*config -> snow_velocity*delta_t
								*grid -> area[i + NO_OF_VECTORS - NO_OF_SCALARS_H]/grid -> volume[base_index];
							}
							// rain
							else if (k < NO_OF_CONDENSED_CONSTITUENTS/2)
							{
								d_vector[j] += impl_weight*config -> rain_velocity*delta_t
								*grid -> area[i + NO_OF_VECTORS - NO_OF_SCALARS_H]/grid -> volume[base_index];
							}
							// clouds
							else if (k < NO_OF_CONDENSED_CONSTITUENTS)
							{
								d_vector[j] += impl_weight*config -> cloud_droplets_velocity*delta_t
								*grid -> area[i + NO_OF_VECTORS - NO_OF_SCALARS_H]/grid -> volume[base_index];
							}
						}
						else
						{
							d_vector[j] = 1.0;
							if (vertical_flux_vector_impl[j - 1] >= 0.0)
							{
								d_vector[j] += impl_weight*delta_t/grid -> volume[base_index]*vertical_flux_vector_impl[j - 1];
							}
							if (vertical_flux_vector_impl[j] < 0.0)
							{
								d_vector[j] -= impl_weight*delta_t/grid -> volume[base_index]*vertical_flux_vector_impl[j];	
							}
						}
						// the explicit component
						// mass densities
						if (quantity_id == 0)
						{
							r_vector[j] =
							state_old -> rho[k*NO_OF_SCALARS + base_index]
							+ delta_t*state_tendency -> rho[k*NO_OF_SCALARS + base_index];
						}
						// density x temperatures
						if (quantity_id == 1)
						{
							r_vector[j] =
							state_old -> condensed_density_temperatures[k*NO_OF_SCALARS + base_index]
							+ delta_t*state_tendency -> condensed_density_temperatures[k*NO_OF_SCALARS + base_index];
						}
						// adding the explicit part of the vertical flux divergence
						if (j == 0)
						{
							r_vector[j] += expl_weight*delta_t*vertical_flux_vector_rhs[j]/grid -> volume[base_index];
						}
						else if (j == NO_OF_LAYERS - 1)
						{
							r_vector[j] += -expl_weight*delta_t*vertical_flux_vector_rhs[j - 1]/grid -> volume[base_index];
							// precipitation
							// snow
							if (k < NO_OF_CONDENSED_CONSTITUENTS/4)
							{
								r_vector[j] += -expl_weight*config -> snow_velocity*delta_t*state_old -> rho[k*NO_OF_SCALARS + i + NO_OF_SCALARS - NO_OF_SCALARS_H]
								*grid -> area[i + NO_OF_VECTORS - NO_OF_SCALARS_H]/grid -> volume[base_index];
							}
							// rain
							else if (k < NO_OF_CONDENSED_CONSTITUENTS/2)
							{
								r_vector[j] += -expl_weight*config -> rain_velocity*delta_t*state_old -> rho[k*NO_OF_SCALARS + i + NO_OF_SCALARS - NO_OF_SCALARS_H]
								*grid -> area[i + NO_OF_VECTORS - NO_OF_SCALARS_H]/grid -> volume[base_index];
							}
							// clouds
							else if (k < NO_OF_CONDENSED_CONSTITUENTS)
							{
								r_vector[j] += -expl_weight*config -> cloud_droplets_velocity*delta_t*state_old -> rho[k*NO_OF_SCALARS + i + NO_OF_SCALARS - NO_OF_SCALARS_H]
								*grid -> area[i + NO_OF_VECTORS - NO_OF_SCALARS_H]/grid -> volume[base_index];
							}
						}
						else
						{
							r_vector[j] += expl_weight*delta_t*(-vertical_flux_vector_rhs[j - 1] + vertical_flux_vector_rhs[j])/grid -> volume[base_index];
						}
					}
					
					// calling the algorithm to solve the system of linear equations
					thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, NO_OF_LAYERS);
					
					// this should account for round-off errors only
					for (int j = 0; j < NO_OF_LAYERS; ++j)
					{
						if (solution_vector[j] < 0.0)
						{
							solution_vector[j] = 0.0;
						}
					}
					
					// writing the result into the new state
					for (int j = 0; j < NO_OF_LAYERS; ++j)
					{
						base_index = i + j*NO_OF_SCALARS_H;
						// mass densities
						if (quantity_id == 0)
						{
							state_new -> rho[k*NO_OF_SCALARS + base_index] = solution_vector[j];
						}
						
						// density x temperature fields
						if (quantity_id == 1)
						{
							state_new -> condensed_density_temperatures[k*NO_OF_SCALARS + base_index] = solution_vector[j];
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
	e_prime_vector[0] = e_vector[0]/d_vector[0];
	for (int j = 1; j < solution_length - 1; ++j)
	{
		e_prime_vector[j] = e_vector[j]/(d_vector[j] - e_prime_vector[j - 1]*c_vector[j - 1]);
	}
	// downward sweep (right-hand side)
	r_prime_vector[0] = r_vector[0]/d_vector[0];
	for (int j = 1; j < solution_length; ++j)
	{
		r_prime_vector[j] = (r_vector[j] - r_prime_vector[j - 1]*c_vector[j - 1])/(d_vector[j] - e_prime_vector[j - 1]*c_vector[j - 1]);
	}
	
	// upward sweep (final solution)
	solution_vector[solution_length - 1] = r_prime_vector[solution_length - 1];
	for (int j = solution_length - 2; j >= 0; --j)
	{
		solution_vector[j] = r_prime_vector[j] - e_prime_vector[j]*solution_vector[j + 1];
	}
	return 0;
}










