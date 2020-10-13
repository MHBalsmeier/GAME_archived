/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../diagnostics/diagnostics.h"

int momentum_diff_diss(State *state, Diagnostics *diagnostics, Diffusion_info *diffusion, Config_info *config_info, Grid *grid)
{	
	// Firstly the diffusion.
	// Evaluating necessary differential operators.
	divv_h(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
	add_vertical_divv(state -> velocity_gas, diagnostics -> velocity_gas_divv, grid);
    grad(diagnostics -> velocity_gas_divv, diffusion -> friction_acc, grid);
    curl_of_vorticity_m(diagnostics -> rel_vort, diagnostics -> curl_of_vorticity_m, grid);
    // Calculating the effective viscosity coefficients.
    calc_divv_term_viscosity_eff(state, config_info, diffusion -> divv_term_viscosity_eff);
    calc_curl_term_viscosity_eff(state, config_info, diffusion -> curl_term_viscosity_eff);
    // Multiplying the values of the differential operators by the effective viscosities.
	scalar_times_vector(diffusion -> divv_term_viscosity_eff, diffusion -> friction_acc, diffusion -> friction_acc, grid, 0);
	scalar_times_vector(diffusion -> curl_term_viscosity_eff, diagnostics -> curl_of_vorticity_m, diffusion -> friction_acc, grid, 1);
	// Then the dissipation (preliminary version).
	// to be implemented
	// Now, the isobaric specific heat capacity of the humid air (the gas phase) needs to be diagnozed.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		if (config_info -> tracers_on == 1)
			diagnostics -> c_h_v_field[i] = spec_heat_cap_diagnostics_v(state -> density_dry[i], state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
		else
			diagnostics -> c_h_v_field[i] = C_D_V;
	}
	 // The heating is a power density, as such it has the SI unit J/(m^3s). As a consequence, this for loop is necessary.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diffusion -> heating_diss[i] = diagnostics -> c_h_v_field[i]*(state -> density_dry[i] + state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i])*diffusion -> heating_diss[i];
	}
	return 0;
}

int curl_of_vorticity_m(Curl_field vorticity, Vector_field out_field, Grid *grid)
{
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_VECTORS; ++i)
	{
		out_field[i] = 0;
	}
	return 0;
}




