/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../integrators.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../../diagnostics/diagnostics.h"

int backward_tendencies(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int no_rk_step)
{
    // Radiation is updated here.
    if (config_info -> rad_update == 1)
    {
    	// radiation will go here
        ;
    }
    
    // Temperature diffusion gets updated here, but only at the last RK step and if heat conduction is switched on.
    if (no_rk_step == 2 && (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1))
    {
    	// Now we need to calculate the scalar diffusion coefficients.
        calc_temp_diffusion_coeffs(state_old, config_info, irreversible_quantities -> scalar_diffusion_coeff_numerical_h, irreversible_quantities -> scalar_diffusion_coeff_numerical_v);
        // The diffusion of the temperature depends on its gradient.
		grad(state_old -> temperature_gas, diagnostics -> temperature_gradient, grid);
		// Now the diffusive temperature flux density can be obtained.
        scalar_times_vector_scalar_h_v(irreversible_quantities -> scalar_diffusion_coeff_numerical_h, irreversible_quantities -> scalar_diffusion_coeff_numerical_v, diagnostics -> temperature_gradient, diagnostics -> flux_density, grid);
        // The divergence of the diffusive temperature flux density is the diffusive temperature heating.
        divv_h(diagnostics -> flux_density, irreversible_quantities -> temperature_diffusion_heating, grid);
    }
    
	integrate_generalized_densities(state_old, state_new, interpolation, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, no_rk_step);
    return 0;
}
