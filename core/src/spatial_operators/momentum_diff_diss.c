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
	divv_h(state -> velocity_gas, diagnostics -> velocity_gas_divv_h, grid);
        grad(diagnostics -> velocity_gas_divv_h, diffusion -> friction_acc, grid);
        // Diagnozing the effective (including turbulence) shear viscosity coefficient.
        calc_kinematic_shear_viscosity_eff(state, config_info, diffusion -> kinematic_shear_viscosity_eff);
    	scalar_times_vector(diffusion -> kinematic_shear_viscosity_eff, diffusion -> friction_acc, diffusion -> friction_acc, grid);
    	
    	// Then the dissipation (preliminary version).
	inner_product(state -> velocity_gas, diffusion -> friction_acc, diffusion -> heating_diss, grid);
	// Now, diagnozing the isobaric specific heat capacity of the humid air (the gas phase) needs to be carried out. 
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
		diffusion -> heating_diss[i] = -diagnostics -> c_h_v_field[i]*(state -> density_dry[i] + state -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i])*diffusion -> heating_diss[i];
	}
	return 0;
}
